"""Stage 3 — joint K × T decision matrix and (K*, T*) recommendation.

Reads:

* Stage-1 per-T outputs (top-5 K=3 and K=4 K-sets, spread metadata,
  Spearman correlations) from
  ``outputs/02_calibration/t_sensitivity_historical/T{T:02d}/``.
* Stage-2 KS / Frobenius / coverage diagnostics from
  ``outputs/02_calibration/t_sensitivity_kirsch_compare/``.
* Stage-3 timing.json files from
  ``outputs/02_calibration/eval_cost_timing/T{T:02d}/``.

Produces:

* ``decision_matrix.csv`` — long-form (K-set, T) → 5-component score
  vector + composite score.
* ``decision_matrix_heatmap_K{K}.csv`` — wide-form table per K (rows =
  ranks 1..5, columns = surviving T values, cell = composite score).
* ``pareto_front_KxT.json`` — non-dominated (K-set, T) tuples on the
  5 normalised criteria + the recommended (K*, T*) and two alternates.
* ``eval_cost_vs_T.csv`` — fit and projection of MOEA wall-time vs T.

Composite score per cell ∈ [0, 1] is the geometric-mean of:

1. ``min_robust_spread`` of the K-set on historical T-blocks at T,
   normalised by the matrix-wide max.
2. ``concept_diversity`` (always 1 if strict rung used, else 1/(rung+1)).
3. ``1 − max_pairwise_|ρ_S|`` of the K-set at T.
4. ``1 − KS_max`` of the K-set's metrics in the Kirsch ensemble at T.
5. ``1 − norm_compute_cost(T)`` from the timing fit.

A geometric mean penalises any low component (don't recommend a K-set
that is fast and low-correlation but has poor Kirsch fidelity).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.extended import (  # noqa: E402
    CONCEPT_MAP,
    HAZARD_CLEAN_METRICS,
)
from src.metrics.screening import (  # noqa: E402
    cluster_metrics as _cluster_metrics,
    enumerate_k_sets,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
HIST_DRIVER = "t_sensitivity_historical"
KIRSCH_COMPARE_DRIVER = "t_sensitivity_kirsch_compare"
TIMING_DRIVER = "eval_cost_timing"
AGG_DRIVER = "t_sensitivity_aggregate"
DRIVER = "decision_matrix"

#: MOEA-FIND production NFE (for cost projection only).
PROD_NFE = 200_000

#: MOEA-FIND production rank count (for cost projection only).
PROD_RANKS = 120

#: Maximum K-set rank to consider per T (top N from Stage-1 alternatives).
MAX_K_RANK = 5

#: Per-T zero-degeneracy gates (excluded before re-enumeration). A metric
#: with frac_zero > GATE_ZERO at a given T is unusable as a Borg axis
#: there because Borg can't move along an axis whose historical
#: reference is zero in too many blocks (e.g. SSI-12 metrics at T=5
#: where ~43% of 5-yr blocks have no SSI-12 events).
DEGENERACY_FRAC_ZERO_GATE = 0.25
DEGENERACY_FRAC_NAN_GATE = 0.10

#: Saturation gate. A metric with > GATE_SAT of its historical T-blocks
#: pinned at the per-T maximum value cannot serve as a useful Borg axis
#: because the historical reference distribution collapses against the
#: upper bound — Borg cannot push the synthetic ensemble beyond a
#: saturated axis without the axis itself becoming non-informative.
#: ``worst_severity`` is a canonical example: 42% sat at T=5, 96% at
#: T=20, 100% at T=30.
DEGENERACY_FRAC_SAT_GATE = 0.25

#: Slug pattern for the Kirsch ensemble produced by Stage 2 / Tier-G
#: recharacterise. Used to load the synthetic-space correlation matrix
#: when ``--correlation-source kirsch`` is selected.
KIRSCH_SLUG_FMT = "n{n_traces}_t{T}_ssi3-12_s{seed}"

#: Concepts to exclude from K-set enumeration by default. Trend metrics
#: (e.g. ``sen_slope_annual_min_neg``) measure a within-block
#: monotonic shift, not a steady-state hazard characteristic — they
#: are not interpretable for stakeholders as "drought hazard" axes
#: even though they produce high robust spread on the historical
#: T-block matrix.
DEFAULT_EXCLUDE_CONCEPTS = ["trend"]

#: Constraint-rung penalty: strict=1, distinct_clusters_only=0.7,
#: distinct_concepts_only=0.7, pairwise_correlation_only=0.4.
RUNG_PENALTY = {
    "distinct_clusters_and_concepts": 1.0,
    "distinct_clusters_only": 0.7,
    "distinct_concepts_only": 0.7,
    "pairwise_correlation_only": 0.4,
    "none": 0.0,
}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _hist_dir(T: int) -> Path:
    return stage_output_dir(STAGE, HIST_DRIVER, slug=f"T{T:02d}", create=False)


def _surviving_T_grid_kirsch() -> Optional[List[int]]:
    p = stage_output_dir(STAGE, KIRSCH_COMPARE_DRIVER, create=False) \
        / "surviving_T_grid_kirsch.json"
    if not p.exists():
        return None
    return list(json.loads(p.read_text()).get("surviving_T_grid", []))


def _surviving_T_grid_hist() -> Optional[List[int]]:
    p = stage_output_dir(STAGE, AGG_DRIVER, create=False) \
        / "surviving_T_grid.json"
    if not p.exists():
        return None
    return list(json.loads(p.read_text()).get("surviving_T_grid", []))


def _load_top_ksets(T: int, K: int) -> Tuple[List[Dict], str]:
    p = _hist_dir(T) / f"k{K}_alternatives.json"
    if not p.exists():
        return [], "none"
    raw = json.loads(p.read_text())
    return raw.get("alternatives", [])[:MAX_K_RANK], raw.get("constraint_rung", "none")


def _load_kirsch_correlations(
    T: int,
    n_traces: int = 10_000,
    seed: int = 42,
    metrics: Optional[List[str]] = None,
) -> Optional[pd.DataFrame]:
    """Load the Stage-2 Kirsch ensemble characteristics and compute the
    Spearman correlation matrix on it (subsetted to ``metrics``).

    Returns ``None`` if the artefact is missing.
    """
    p = stage_output_dir(
        "03_kirsch_library", "build_library_extended",
        slug=KIRSCH_SLUG_FMT.format(n_traces=n_traces, T=T, seed=seed),
        create=False,
    ) / "characteristics_extended.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    df = pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])
    if metrics is not None:
        cols = [m for m in metrics if m in df.columns]
        df = df[cols]
    return df.astype(float).corr(method="spearman")


def _load_per_T_screening_artifacts(T: int):
    """Return (spread_df, cluster_df, spearman_corr, surviving_metrics)
    for one T value, ready to feed
    :func:`src.metrics.screening.enumerate_k_sets` directly. Applies the
    per-T zero-degeneracy filter on top of the existing screen-pass
    flag so SSI-12 metrics that are zero in a large fraction of
    historical blocks are excluded BEFORE enumeration (rather than
    filtering them out post-hoc).
    """
    d = _hist_dir(T)
    spread = pd.read_csv(d / "per_metric_spread.csv")
    deg = pd.read_csv(d / "degeneracy.csv")
    spearman = pd.read_csv(d / "spearman_corr.csv", index_col=0)
    cluster_df = pd.read_json(d / "clusters.json")

    deg_lookup = deg.set_index("metric")
    survivors: List[str] = []
    for _, row in spread.iterrows():
        if not bool(row["passes_screen"]):
            continue
        m = row["metric"]
        if m in deg_lookup.index:
            row_deg = deg_lookup.loc[m]
            frac_zero = float(row_deg.get("frac_zero") or 0.0)
            frac_nan = float(row_deg.get("frac_nan") or 0.0)
            frac_sat = float(row_deg.get("frac_saturated_at_max") or 0.0)
            if frac_zero > DEGENERACY_FRAC_ZERO_GATE:
                continue
            if frac_nan > DEGENERACY_FRAC_NAN_GATE:
                continue
            if frac_sat > DEGENERACY_FRAC_SAT_GATE:
                continue
        if m in spearman.columns:
            survivors.append(m)
    return spread, cluster_df, spearman, survivors


def _enumerate_filtered_ksets(
    T: int,
    K: int,
    correlation_cap: float,
    exclude_concepts: List[str],
    n_alternatives: int,
    *,
    metric_allowlist: Optional[List[str]] = None,
    correlation_source: str = "kirsch",
    n_traces: int = 10_000,
    seed: int = 42,
) -> Tuple[List[Dict], str]:
    """Re-enumerate K-sets at this T with hazard-only filter + per-T
    zero-degeneracy filter. Falls back to the Stage-1 enumeration
    artefact if any input is missing.
    """
    try:
        spread, cluster_df, spearman, surviving = \
            _load_per_T_screening_artifacts(T)
    except FileNotFoundError:
        return _load_top_ksets(T, K)

    if exclude_concepts:
        excl = set(exclude_concepts)
        cluster_df = cluster_df[~cluster_df["concept"].isin(excl)]
        surviving = [m for m in surviving
                     if CONCEPT_MAP.get(m, "?") not in excl]
        spearman = spearman.loc[surviving, surviving] if surviving else spearman

    if metric_allowlist is not None:
        allowed = set(metric_allowlist)
        cluster_df = cluster_df[cluster_df["metric"].isin(allowed)]
        surviving = [m for m in surviving if m in allowed]
        spearman = spearman.loc[surviving, surviving] if surviving else spearman

    if len(surviving) < K:
        return [], "none"

    # Optionally swap the historical correlation matrix for the Kirsch
    # ensemble's. The Kirsch matrix has 10k independent observations
    # (vs the 14 effective independent samples from a 73-WY record at
    # T=5), so it's the statistically tighter and operationally
    # correct estimator for inter-metric correlations Borg will see
    # while searching over the synthetic ensemble.
    if correlation_source in ("kirsch", "both") and surviving:
        kirsch_corr = _load_kirsch_correlations(
            T, n_traces=n_traces, seed=seed, metrics=surviving,
        )
        if kirsch_corr is None:
            print(f"[warn] T={T}: Kirsch correlations unavailable; "
                  f"falling back to historical")
        else:
            kirsch_corr = kirsch_corr.reindex(
                index=surviving, columns=surviving,
            )
            if correlation_source == "kirsch":
                # Recluster on the Kirsch correlation skeleton so cluster
                # membership reflects synthetic-space redundancy.
                spread_lookup = spread.set_index("metric")["spread_score"]
                cluster_df = _cluster_metrics(kirsch_corr, spread_lookup)
                if metric_allowlist is not None:
                    cluster_df = cluster_df[
                        cluster_df["metric"].isin(set(metric_allowlist))
                    ]
                spearman = kirsch_corr
            else:  # "both" — use the elementwise max of |hist| and |kirsch|
                hist_abs = spearman.abs()
                kir_abs = kirsch_corr.abs()
                # Reuse the union for the cap filter, keeping signs from hist
                signs = np.sign(spearman.values)
                merged = np.maximum(hist_abs.values, kir_abs.values)
                spearman = pd.DataFrame(
                    merged * signs,
                    index=spearman.index, columns=spearman.columns,
                )

    rungs = [
        ("distinct_clusters_and_concepts", True, True),
        ("distinct_clusters_only", True, False),
        ("distinct_concepts_only", False, True),
        ("pairwise_correlation_only", False, False),
    ]
    for label, dc, dcp in rungs:
        cands = enumerate_k_sets(
            surviving, spearman, cluster_df, spread, K,
            correlation_cap=correlation_cap,
            require_distinct_clusters=dc,
            require_distinct_concepts=dcp,
        )
        if cands:
            return cands[:n_alternatives], label
    return [], "none"


def _load_ks_table() -> pd.DataFrame:
    p = stage_output_dir(STAGE, KIRSCH_COMPARE_DRIVER, create=False) \
        / "ks_vs_historical.csv"
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p)


def _load_timing(T: int) -> Optional[Dict]:
    p = stage_output_dir(STAGE, TIMING_DRIVER, slug=f"T{T:02d}", create=False) \
        / "timing.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


# ---------------------------------------------------------------------------
# Cost model
# ---------------------------------------------------------------------------


def fit_cost_model(timings: Dict[int, Dict]) -> Dict[str, float]:
    """Linear fit of median per-evaluation wall-time vs T."""
    if len(timings) < 2:
        return {"a": float("nan"), "b": float("nan"),
                "T_grid": list(timings.keys()),
                "median_eval_seconds": [t["median"] for t in timings.values()]}
    Ts = np.array(sorted(timings.keys()), dtype=float)
    medians = np.array([timings[int(t)]["median"] for t in Ts])
    A = np.vstack([np.ones_like(Ts), Ts]).T
    sol, *_ = np.linalg.lstsq(A, medians, rcond=None)
    a, b = float(sol[0]), float(sol[1])
    pred = a + b * Ts
    ss_res = float(np.sum((medians - pred) ** 2))
    ss_tot = float(np.sum((medians - medians.mean()) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    return {
        "a_seconds": a,
        "b_seconds_per_year": b,
        "r2": r2,
        "T_grid": [int(t) for t in Ts],
        "median_eval_seconds": medians.tolist(),
    }


def projected_walltime_hours(T: int, fit: Dict[str, float]) -> float:
    """Project total MOEA wall-time at T using the cost-model fit."""
    a = fit.get("a_seconds")
    b = fit.get("b_seconds_per_year")
    if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
        return float("nan")
    eval_s = a + b * T
    total_s = eval_s * PROD_NFE / PROD_RANKS
    return float(total_s / 3600.0)


# ---------------------------------------------------------------------------
# Decision matrix
# ---------------------------------------------------------------------------


def build_decision_matrix(
    T_grid: List[int],
    ks_table: pd.DataFrame,
    fit: Dict[str, float],
    *,
    correlation_cap: float,
    exclude_concepts: List[str],
    n_alternatives: int,
    metric_allowlist: Optional[List[str]] = None,
    correlation_source: str = "kirsch",
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    cost_hours = {T: projected_walltime_hours(T, fit) for T in T_grid}

    for T in T_grid:
        for K in (3, 4):
            alts, rung = _enumerate_filtered_ksets(
                T, K,
                correlation_cap=correlation_cap,
                exclude_concepts=exclude_concepts,
                n_alternatives=n_alternatives,
                metric_allowlist=metric_allowlist,
                correlation_source=correlation_source,
            )
            for rank, alt in enumerate(alts, 1):
                metrics = list(alt["metrics"])
                ks_subset = ks_table[
                    (ks_table["T_years"] == T)
                    & (ks_table["metric"].isin(metrics))
                ]
                ks_max = float(ks_subset["ks_statistic"].max()) \
                    if not ks_subset.empty else float("nan")
                rows.append({
                    "T_years": T,
                    "K": K,
                    "rank": rank,
                    "metrics": "|".join(metrics),
                    "concepts": "|".join(alt["concepts"]),
                    "constraint_rung": rung,
                    "rung_penalty": RUNG_PENALTY.get(rung, 0.0),
                    "min_spread": float(alt["min_spread"]),
                    "sum_spread": float(alt["sum_spread"]),
                    "max_pairwise_rho": float(alt["max_pairwise_rho"]),
                    "ks_max_metric": ks_max,
                    "projected_walltime_hours_120rank": cost_hours[T],
                })
    return pd.DataFrame(rows)


def normalise_and_score(df: pd.DataFrame) -> pd.DataFrame:
    """Add normalised components and the composite score.

    Components are clipped to [eps, 1] so the geometric mean is finite
    even when one component is exactly zero.
    """
    if df.empty:
        return df
    out = df.copy()

    # 1) min spread (higher is better)
    s = out["min_spread"].astype(float)
    s_max = s.max()
    out["score_spread"] = (s / s_max).clip(lower=1e-3) \
        if s_max > 0 else np.nan

    # 2) constraint-rung penalty (already in [0, 1])
    out["score_rung"] = out["rung_penalty"].astype(float).clip(lower=1e-3)

    # 3) low correlation: 1 − max |ρ|
    out["score_indep"] = (1.0 - out["max_pairwise_rho"]).clip(lower=1e-3)

    # 4) Kirsch fidelity: 1 − KS_max
    ks = out["ks_max_metric"].astype(float)
    out["score_fidelity"] = (1.0 - ks.fillna(1.0)).clip(lower=1e-3)

    # 5) cost: 1 − cost / max_cost
    c = out["projected_walltime_hours_120rank"].astype(float)
    c_max = c.max()
    out["score_cost"] = (1.0 - c / c_max).clip(lower=1e-3) \
        if c_max > 0 and np.isfinite(c_max) else 1.0

    components = ["score_spread", "score_rung", "score_indep",
                  "score_fidelity", "score_cost"]
    out["composite_score"] = np.exp(
        out[components].apply(np.log).mean(axis=1)
    )
    return out


def pareto_front(df: pd.DataFrame) -> pd.DataFrame:
    """Non-dominated (K-set, T) rows on the 5 score components.

    Maximisation on every score_* column. O(N²) — fine for ≤40 rows.
    """
    cols = ["score_spread", "score_rung", "score_indep",
            "score_fidelity", "score_cost"]
    if df.empty:
        return df
    M = df[cols].values
    n = M.shape[0]
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        for j in range(n):
            if i == j or not keep[j]:
                continue
            if np.all(M[j] >= M[i]) and np.any(M[j] > M[i]):
                keep[i] = False
                break
    return df[keep].copy().sort_values("composite_score", ascending=False)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-grid", type=int, nargs="+", default=None,
                   help="Override the T grid. Default uses the Stage-1 "
                        "hist survivor list (Kirsch fidelity is informational).")
    p.add_argument("--correlation-cap", type=float, default=0.6,
                   help="Pairwise |Spearman ρ| cap for K-set enumeration.")
    p.add_argument("--exclude-concepts", type=str, nargs="*",
                   default=DEFAULT_EXCLUDE_CONCEPTS,
                   help="Hydrologic concepts to exclude from K-set "
                        "enumeration (default: trend). Pass empty list "
                        "with `--exclude-concepts ` to include all "
                        "concepts including within-block trend.")
    p.add_argument("--n-alternatives", type=int, default=10,
                   help="Top-N K-sets to surface per (T, K). Higher than "
                        "Stage-1's 5 so the matrix has more candidates.")
    p.add_argument("--metric-pool", choices=["full", "hazard-clean"],
                   default="hazard-clean",
                   help="Metric pool fed into K-set enumeration. "
                        "'hazard-clean' restricts to the Tier-G + "
                        "time_in_drought_fraction set per DD-15 "
                        "(continuous, full-record-anchored, non-averaged, "
                        "non-stratified). 'full' uses all 34 candidates "
                        "subject to the saturation/zero/concept filters.")
    p.add_argument("--correlation-source",
                   choices=["hist", "kirsch", "both"],
                   default="kirsch",
                   help="Spearman correlation reference for K-set "
                        "redundancy filtering. 'kirsch' (default, DD-15a) "
                        "uses the 10k-trace synthetic ensemble — the "
                        "correlation Borg actually sees. 'hist' uses the "
                        "stride-1 historical blocks (heavily inflated by "
                        "overlap at small T). 'both' takes the elementwise "
                        "max of |ρ_hist| and |ρ_kirsch| for the cap filter.")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/{DRIVER}] → {out_dir}")

    # --- Determine surviving T grid ---
    # Use the historical-block viability list as the primary gate
    # (whether Stage 1 produced strict-rung K-sets at each T). The
    # Stage-2 Kirsch fidelity gate is informational only: it feeds
    # the composite score's `score_fidelity` component below, but
    # does NOT filter T values, because the historical T-block
    # distribution is degenerate by construction at large T (blocks
    # overlap heavily with N_hist = 73 → at T = 30 most blocks share
    # the record-worst drought event, so e.g. `worst_severity`
    # saturates and a 10k-realization Kirsch ensemble sampling that
    # axis legitimately appears "different" by KS without that being
    # a generator deficiency).
    hist_grid = _surviving_T_grid_hist()
    kirsch_grid = _surviving_T_grid_kirsch()
    if args.T_grid is not None:
        T_grid = sorted(set(args.T_grid))
        print(f"[diag] T grid (overridden): {T_grid}")
    elif hist_grid is not None:
        T_grid = sorted(hist_grid)
        kirsch_str = (f"{kirsch_grid}" if kirsch_grid is not None
                      else "(absent)")
        print(f"[diag] hist survivors={hist_grid} (gating); "
              f"kirsch survivors={kirsch_str} (informational)")
    else:
        T_grid = [5, 10, 20, 30]
        print(f"[diag] no survivor lists; defaulting to {T_grid}")

    if not T_grid:
        raise SystemExit("Empty T grid from Stage 1 viability gate.")

    # --- Cost model ---
    timings: Dict[int, Dict] = {}
    for T in T_grid:
        t = _load_timing(T)
        if t is not None:
            timings[T] = t
        else:
            print(f"[warn] T={T}: missing timing.json; cost component will be NaN.")
    fit = fit_cost_model(timings)
    cost_rows = [{
        "T_years": T,
        "median_eval_seconds": timings[T]["median"] if T in timings else float("nan"),
        "mean_eval_seconds": timings[T]["mean"] if T in timings else float("nan"),
        "projected_walltime_hours_120rank": projected_walltime_hours(T, fit),
    } for T in T_grid]
    pd.DataFrame(cost_rows).to_csv(out_dir / "eval_cost_vs_T.csv", index=False)
    (out_dir / "cost_model_fit.json").write_text(json.dumps(fit, indent=2))
    if np.isfinite(fit.get("a_seconds", np.nan)):
        print(f"[diag] cost fit: t_eval(T) = {fit['a_seconds']:.4f} + "
              f"{fit['b_seconds_per_year']:.4f}·T  "
              f"(R²={fit['r2']:.3f})")

    # --- Decision matrix ---
    ks_table = _load_ks_table()
    metric_allowlist = (list(HAZARD_CLEAN_METRICS)
                        if args.metric_pool == "hazard-clean" else None)
    print(f"[diag] re-enumerating K-sets with "
          f"metric_pool={args.metric_pool}, "
          f"correlation_source={args.correlation_source}, "
          f"exclude_concepts={args.exclude_concepts}, "
          f"correlation_cap={args.correlation_cap}, "
          f"n_alternatives={args.n_alternatives}")
    if metric_allowlist is not None:
        print(f"[diag] hazard-clean allowlist ({len(metric_allowlist)}): "
              f"{metric_allowlist}")
    df = build_decision_matrix(
        T_grid, ks_table, fit,
        correlation_cap=args.correlation_cap,
        exclude_concepts=list(args.exclude_concepts),
        n_alternatives=args.n_alternatives,
        metric_allowlist=metric_allowlist,
        correlation_source=args.correlation_source,
    )
    if df.empty:
        raise SystemExit("Decision matrix empty — check Stage-1 K-set artifacts.")
    df = normalise_and_score(df)
    df.to_csv(out_dir / "decision_matrix.csv", index=False)
    print(f"[diag] wrote decision_matrix.csv ({len(df)} rows)")

    for K in (3, 4):
        sub = df[df["K"] == K]
        if sub.empty:
            continue
        wide = sub.pivot_table(
            index="rank", columns="T_years", values="composite_score",
            aggfunc="first",
        )
        wide.to_csv(out_dir / f"decision_matrix_heatmap_K{K}.csv")

    pareto_df = pareto_front(df)
    pareto_df.to_csv(out_dir / "pareto_front_KxT.csv", index=False)

    # --- Recommendation: top by composite, with two alternates ---
    df_sorted = df.sort_values("composite_score", ascending=False)
    top = df_sorted.iloc[0]
    alternates = df_sorted.iloc[1:].copy()
    # Prefer alternates with different (K, T) than the top pick.
    alternates = alternates[
        (alternates["K"] != top["K"]) | (alternates["T_years"] != top["T_years"])
    ].head(2)
    rec = {
        "recommended": {
            "K": int(top["K"]),
            "T_years": int(top["T_years"]),
            "rank_within_T": int(top["rank"]),
            "metrics": top["metrics"].split("|"),
            "concepts": top["concepts"].split("|"),
            "constraint_rung": top["constraint_rung"],
            "composite_score": float(top["composite_score"]),
            "score_spread": float(top["score_spread"]),
            "score_rung": float(top["score_rung"]),
            "score_indep": float(top["score_indep"]),
            "score_fidelity": float(top["score_fidelity"]),
            "score_cost": float(top["score_cost"]),
            "min_spread": float(top["min_spread"]),
            "max_pairwise_rho": float(top["max_pairwise_rho"]),
            "ks_max_metric": float(top["ks_max_metric"])
                if pd.notna(top["ks_max_metric"]) else None,
            "projected_walltime_hours_120rank":
                float(top["projected_walltime_hours_120rank"])
                if pd.notna(top["projected_walltime_hours_120rank"]) else None,
        },
        "alternates": [{
            "K": int(r["K"]),
            "T_years": int(r["T_years"]),
            "metrics": r["metrics"].split("|"),
            "composite_score": float(r["composite_score"]),
        } for _, r in alternates.iterrows()],
        "n_pareto_front": int(len(pareto_df)),
        "T_grid_evaluated": T_grid,
        "cost_model": fit,
    }
    (out_dir / "pareto_front_KxT.json").write_text(json.dumps(rec, indent=2))

    print(f"[rec] recommended: K={top['K']}, T={top['T_years']}, "
          f"metrics={top['metrics']}, score={top['composite_score']:.3f}")
    print("[diag] done.")


if __name__ == "__main__":
    main()
