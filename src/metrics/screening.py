"""Shared drought-metric screening API.

Extracted from ``workflows/02_calibration/metric_explorer.py`` so that
both the original single-T driver and the new T-sweep drivers in
``workflows/02_calibration/t_sensitivity_*.py`` share a single
implementation. Pure relocation — no logic change.

Public surface:

* Constants: :data:`CLUSTER_DISTANCE_CUT`, :data:`SKEW_DROP_THRESHOLD`,
  :data:`MIN_STD`, :data:`N_ALTERNATIVES`.
* :func:`compute_block_matrix` — run all 28 candidate metrics on every
  T-year block.
* :func:`per_metric_spread` — IQR / range / robust spread / skew per
  metric, with screen-pass flag.
* :func:`correlation_matrices` — Pearson and Spearman matrices on a
  surviving subset.
* :func:`cluster_metrics` — average-linkage hierarchical clustering on
  ``1 − |ρ_S|``, with representative chosen by max robust spread.
* :func:`enumerate_k_sets` — brute-force K-subset enumeration under the
  strict clusters-AND-concepts constraints + pairwise |ρ_S| cap.
* :func:`relax_until_nonempty` — strict-then-relaxed enumeration ladder.
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from src.metrics.extended import (
    CANDIDATE_METRIC_NAMES,
    CONCEPT_MAP,
    FullRecordRefs,
    compute_all_candidates,
)
from src.hydrology.historical_blocks import resample_historical_blocks

#: Cut threshold on ``1 − |Spearman ρ|`` distance (i.e. metrics with
#: |ρ| ≥ 0.7 are clustered as redundant).
CLUSTER_DISTANCE_CUT = 0.30

#: Skewness magnitude above which a metric is considered too one-sided
#: to serve as a usable MOEA axis.
SKEW_DROP_THRESHOLD = 3.0

#: Minimum standard deviation under which a metric is treated as
#: essentially constant across T-blocks and dropped.
MIN_STD = 1e-9

#: How many top-scoring K-sets to surface in the recommendation report.
N_ALTERNATIVES = 5


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------


def compute_block_matrix(
    monthly_1d: np.ndarray,
    T_years: int,
    stride: int,
    ssi3_calc,
    ssi12_calc,
    q80_threshold: float,
    *,
    full_record_refs: "FullRecordRefs | None" = None,
) -> pd.DataFrame:
    """Run every candidate metric on every T-year block.

    ``full_record_refs`` carries the reference statistics needed by the
    Tier-G hazard-clean metrics (DD-15). Pass ``None`` to fall back to
    the legacy 28-metric pool with Tier G zero-filled.
    """
    blocks_1d = resample_historical_blocks(monthly_1d, T_years, stride)
    rows: List[Dict[str, float]] = []
    for blk in blocks_1d:
        blk_2d = blk.reshape(T_years, 12)
        chars = compute_all_candidates(
            blk, blk_2d, ssi3_calc, ssi12_calc, q80_threshold,
            full_record_refs=full_record_refs,
        )
        rows.append(chars)
    df = pd.DataFrame(rows)
    diag_cols = [c for c in df.columns if c not in CANDIDATE_METRIC_NAMES]
    cand_cols = [c for c in CANDIDATE_METRIC_NAMES if c in df.columns]
    return df[cand_cols + diag_cols]


# ---------------------------------------------------------------------------
# Spread screening
# ---------------------------------------------------------------------------


def per_metric_spread(chars_df: pd.DataFrame) -> pd.DataFrame:
    """Per-metric spread descriptors with robust score and screen-pass flag.

    ``spread_score = IQR / (|median| + σ)`` — stable when ``|median|`` is
    small (sign-crossing metrics like ``sen_slope_annual_min_neg`` no
    longer get artificially inflated rankings). Falls back to ``IQR / σ``
    when the median is exactly zero, and to ``0`` when both are zero.

    Screening drops metrics with std < ``MIN_STD`` or |skew| >
    ``SKEW_DROP_THRESHOLD``.
    """
    rows = []
    for name in CANDIDATE_METRIC_NAMES:
        if name not in chars_df.columns:
            continue
        v = chars_df[name].astype(float).values
        if v.size == 0:
            continue
        mean = float(np.mean(v))
        med = float(np.median(v))
        std = float(np.std(v, ddof=1)) if v.size > 1 else 0.0
        q1, q3 = np.percentile(v, [25.0, 75.0])
        iqr = float(q3 - q1)
        denom = abs(med) + std
        spread_score = float(iqr / denom) if denom > 1e-12 else 0.0
        skew = float(stats.skew(v, bias=False)) if v.size > 2 else 0.0
        # A metric where >50% of blocks share the same value (e.g.
        # ``max_duration`` saturating at 9 months across most 20-yr
        # blocks) is degenerate — Borg can't move along an axis that
        # has only one or two distinct values.
        passes = (
            std >= MIN_STD
            and iqr > 0.0
            and spread_score > 0.0
            and np.isfinite(spread_score)
            and abs(skew) <= SKEW_DROP_THRESHOLD
        )
        rows.append({
            "metric": name,
            "concept": CONCEPT_MAP.get(name, "unknown"),
            "mean": mean,
            "median": med,
            "std": std,
            "min": float(np.min(v)),
            "max": float(np.max(v)),
            "range": float(np.max(v) - np.min(v)),
            "iqr": iqr,
            "spread_score": spread_score,
            "skew": skew,
            "passes_screen": bool(passes),
        })
    return (
        pd.DataFrame(rows)
        .sort_values("spread_score", ascending=False, na_position="last")
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Correlation + hierarchical clustering
# ---------------------------------------------------------------------------


def correlation_matrices(
    chars_df: pd.DataFrame,
    metrics: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Pearson and Spearman correlation matrices for the surviving metrics."""
    sub = chars_df[list(metrics)].astype(float)
    pearson = sub.corr(method="pearson")
    spearman = sub.corr(method="spearman")
    return pearson, spearman


def cluster_metrics(
    spearman_corr: pd.DataFrame,
    spread_score: pd.Series,
    distance_cut: float = CLUSTER_DISTANCE_CUT,
) -> pd.DataFrame:
    """Hierarchical-cluster metrics on ``1 − |ρ_S|`` and pick a representative.

    Average-linkage clustering, cut at ``distance_cut`` (default 0.30 →
    metrics with |ρ_S| ≥ 0.7 share a cluster). Within each cluster the
    representative is the metric with the highest **robust spread score**
    (``IQR / (|median| + σ)``).
    """
    metrics = list(spearman_corr.columns)
    if len(metrics) < 2:
        return pd.DataFrame(
            {"metric": metrics,
             "cluster_id": [1] * len(metrics),
             "concept": [CONCEPT_MAP.get(m, "unknown") for m in metrics],
             "is_representative": [True] * len(metrics)}
        )
    dist = 1.0 - spearman_corr.abs().values
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0.0, None)
    dist = (dist + dist.T) / 2.0
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    cluster_ids = fcluster(Z, t=distance_cut, criterion="distance")

    rows = pd.DataFrame({
        "metric": metrics,
        "cluster_id": cluster_ids.astype(int),
        "concept": [CONCEPT_MAP.get(m, "unknown") for m in metrics],
    })
    rows["spread_score"] = rows["metric"].map(spread_score).fillna(-np.inf)
    rows["is_representative"] = False
    for cid, group in rows.groupby("cluster_id", sort=False):
        winner = group["spread_score"].idxmax()
        rows.loc[winner, "is_representative"] = True
    return rows[["metric", "cluster_id", "concept", "is_representative"]]


# ---------------------------------------------------------------------------
# Brute-force K-set enumeration
# ---------------------------------------------------------------------------


def enumerate_k_sets(
    pool: List[str],
    spearman_corr: pd.DataFrame,
    cluster_df: pd.DataFrame,
    spread_df: pd.DataFrame,
    K: int,
    *,
    correlation_cap: float = 0.6,
    require_distinct_clusters: bool = True,
    require_distinct_concepts: bool = True,
) -> List[Dict]:
    """Enumerate all K-metric subsets of ``pool`` satisfying the hard
    constraints, scored by ``min(spread_score over K members)``.

    Hard constraints:

    * Pairwise |ρ_S| < ``correlation_cap``.
    * No two members share a redundancy cluster (if
      ``require_distinct_clusters``).
    * No two members share a concept tag (if ``require_distinct_concepts``).

    Score: **min spread score** (rewards K-sets where every axis is
    informative — no weak link). Tiebreaker: sum of spread scores.

    Returns the full list of feasible K-sets, sorted by score descending.
    The caller picks the top ``N_ALTERNATIVES`` for the report.
    """
    cluster_lookup = cluster_df.set_index("metric")["cluster_id"].to_dict()
    concept_lookup = cluster_df.set_index("metric")["concept"].to_dict()
    spread_lookup = spread_df.set_index("metric")["spread_score"].to_dict()

    candidates: List[Dict] = []
    for combo in combinations(pool, K):
        if require_distinct_clusters:
            if len({cluster_lookup.get(m) for m in combo}) < K:
                continue
        if require_distinct_concepts:
            if len({concept_lookup.get(m) for m in combo}) < K:
                continue
        sub = spearman_corr.loc[list(combo), list(combo)].values.copy()
        np.fill_diagonal(sub, 0.0)
        max_offdiag = float(np.abs(sub).max())
        if max_offdiag >= correlation_cap:
            continue
        spreads = [float(spread_lookup.get(m, 0.0)) for m in combo]
        candidates.append({
            "metrics": list(combo),
            "concepts": [concept_lookup.get(m, "?") for m in combo],
            "clusters": [int(cluster_lookup.get(m, -1)) for m in combo],
            "spreads": spreads,
            "min_spread": float(min(spreads)),
            "sum_spread": float(sum(spreads)),
            "max_pairwise_rho": max_offdiag,
        })
    candidates.sort(
        key=lambda c: (-c["min_spread"], -c["sum_spread"]),
    )
    return candidates


def relax_until_nonempty(
    pool: List[str],
    spearman_corr: pd.DataFrame,
    cluster_df: pd.DataFrame,
    spread_df: pd.DataFrame,
    K: int,
    *,
    correlation_cap: float,
) -> Tuple[List[Dict], str]:
    """Try the strict-then-relaxed enumeration ladder.

    Tries (1) distinct clusters AND distinct concepts; if no candidates,
    falls back to (2) distinct clusters only; then to (3) distinct concepts
    only. Records which rung produced the result so the report can flag
    when constraints had to be relaxed.
    """
    rungs = [
        ("distinct_clusters_and_concepts", True, True),
        ("distinct_clusters_only", True, False),
        ("distinct_concepts_only", False, True),
        ("pairwise_correlation_only", False, False),
    ]
    for label, dc, dcp in rungs:
        cands = enumerate_k_sets(
            pool, spearman_corr, cluster_df, spread_df, K,
            correlation_cap=correlation_cap,
            require_distinct_clusters=dc,
            require_distinct_concepts=dcp,
        )
        if cands:
            return cands, label
    return [], "none"
