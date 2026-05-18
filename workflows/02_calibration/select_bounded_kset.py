"""DD-15c — diagnostic-driven (mapping × K=4 windows) selection.

Phase 3 of the bounded T=1 reformulation. Consumes the parquets emitted
by ``recompute_bounded_candidates`` and selects the K=4 production
metric set from the 48-candidate pool using objective scoring.

**Per-candidate scores** (combining historical + Pareto archives):

* ``discrimination_hist`` = ``IQR(D over historical T=1) / D*`` — penalises
  windows whose D values collapse to a narrow band.
* ``tail_resolution`` = ``std(D over the top-decile of Pareto by raw s_j)``
  — penalises mappings/windows that flatten the extreme drought tail
  (the DD-15c no-clip requirement). Higher = more tail discrimination.
* ``flood_corner_frac`` = fraction of Pareto members with ``D > 0.9``
  whose ``max_monthly_flow > 10k cfs`` — quantifies the DD-15c flood-
  corner pathology that the new metrics must suppress.
* ``pareto_spread_iqr`` = ``IQR(D over Pareto) / D*`` — avoids degenerate
  "all saturated" candidates.
* ``hist_pareto_w1`` = 1-Wasserstein distance between historical and
  Pareto ECDFs of ``D``; helps rank windows whose Pareto extends the
  drought tail without extending the flood tail.

**Per-K-set scores** (4-element subsets of survivors):

* ``median_abs_rho`` = median pairwise |Spearman ρ| over historical T=1
  blocks. K-sets above 0.6 dropped (existing convention).
* ``concept_diversity`` = number of distinct concept tags in the K-set
  (must be ≥ 4 for full diversity).
* ``composite_score`` = ``√(mean discrimination_hist × (1 − median_abs_rho))``
  — the headline rank.

**Selection algorithm:**

1. Per-candidate gates: ``discrimination_hist >= 0.2``,
   ``tail_resolution >= 1e-4``, ``flood_corner_frac <= 0.05``.
2. Enumerate K=4 subsets per mapping (g vs e) over survivors;
   filter ``median_abs_rho <= 0.6`` and ``concept_diversity >= 4``.
3. Rank by ``composite_score``. Pick best per mapping; compare and
   choose overall winner.

Outputs (under ``outputs/02_calibration/select_bounded_kset/``):

* ``per_candidate_scores.csv`` — table of every candidate × score.
* ``kset_rankings_g.csv`` / ``kset_rankings_e.csv`` — top-N K-sets per
  mapping with composite_score, median_abs_rho, concept tags.
* ``selection_report.json`` — winning (mapping, K=4 windows) tuple plus
  full diagnostic context.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    CANDIDATE_BOUNDED_CONCEPT_MAP,
    CANDIDATE_BOUNDED_METRIC_NAMES,
    WINDOW_SPECS,
)


_FLOOD_THRESHOLD_CFS = 10_000.0
_DSTAR = 1.0  # CONSTANT anti-ideal for all bounded candidates


def _load_recompute_outputs(in_dir: Path) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    hist_path = in_dir / "historical.parquet"
    if not hist_path.exists():
        raise FileNotFoundError(
            f"Phase 2 historical.parquet not found at {hist_path}; "
            f"run recompute_bounded_candidates first."
        )
    hist_df = pd.read_parquet(hist_path)
    archive_dfs: Dict[str, pd.DataFrame] = {}
    for path in sorted(in_dir.glob("*.parquet")):
        if path.name == "historical.parquet":
            continue
        slug = path.stem
        archive_dfs[slug] = pd.read_parquet(path)
    if not archive_dfs:
        raise FileNotFoundError(f"No archive parquets in {in_dir}")
    return hist_df, archive_dfs


def _per_candidate_scores(
    hist_df: pd.DataFrame,
    archive_dfs: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Per-candidate diagnostic scores aggregated across archives."""
    rows: List[Dict[str, float]] = []
    pareto_concat = pd.concat(list(archive_dfs.values()), ignore_index=True)
    for cand in CANDIDATE_BOUNDED_METRIC_NAMES:
        hist_vals = hist_df[cand].to_numpy(dtype=float)
        pareto_vals = pareto_concat[cand].to_numpy(dtype=float)
        # discrimination_hist: IQR / D*
        q25h, q75h = np.percentile(hist_vals, [25, 75])
        discrim = float((q75h - q25h) / _DSTAR)
        # pareto_spread_iqr
        q25p, q75p = np.percentile(pareto_vals, [25, 75])
        spread = float((q75p - q25p) / _DSTAR)
        # tail_resolution: std of D in the top decile of Pareto by D itself
        thresh = float(np.percentile(pareto_vals, 90.0))
        top_decile = pareto_vals[pareto_vals >= thresh]
        tail_res = float(top_decile.std(ddof=0)) if top_decile.size > 1 else 0.0
        # flood_corner_frac: fraction of Pareto with D > 0.9 AND max_monthly_flow > threshold
        high_d_mask = pareto_vals > 0.9
        if high_d_mask.any():
            sub = pareto_concat.loc[high_d_mask, "max_monthly_flow"].to_numpy(dtype=float)
            flood_frac = float((sub > _FLOOD_THRESHOLD_CFS).mean())
        else:
            flood_frac = 0.0
        # 1-Wasserstein distance hist vs pareto
        w1 = float(sp_stats.wasserstein_distance(hist_vals, pareto_vals))
        rows.append({
            "candidate": cand,
            "concept": CANDIDATE_BOUNDED_CONCEPT_MAP.get(cand, "other"),
            "mapping": cand.rsplit("_", 1)[-1],  # last token: 'g' or 'e'
            "window": cand.rsplit("_", 1)[0],
            "discrimination_hist": discrim,
            "pareto_spread_iqr": spread,
            "tail_resolution": tail_res,
            "flood_corner_frac": flood_frac,
            "hist_pareto_w1": w1,
            "n_hist": int(hist_vals.size),
            "n_pareto": int(pareto_vals.size),
        })
    return pd.DataFrame(rows)


def _enumerate_ksets(
    candidates: Iterable[str],
    hist_df: pd.DataFrame,
    *,
    K: int = 4,
) -> pd.DataFrame:
    """Enumerate K-element subsets and score each."""
    cands = list(candidates)
    if len(cands) < K:
        return pd.DataFrame(columns=[
            "metrics", "concepts", "median_abs_rho", "max_abs_rho",
            "n_distinct_concepts", "mean_discrimination", "composite_score",
        ])
    rho_full = hist_df[cands].corr(method="spearman").abs().to_numpy()
    cand_idx = {c: i for i, c in enumerate(cands)}
    rows: List[Dict[str, object]] = []
    for combo in combinations(cands, K):
        idx = [cand_idx[c] for c in combo]
        # Off-diagonal pairwise |rho| within the subset
        sub = rho_full[np.ix_(idx, idx)]
        upper = sub[np.triu_indices_from(sub, k=1)]
        median_rho = float(np.median(upper))
        max_rho = float(np.max(upper))
        concepts = [CANDIDATE_BOUNDED_CONCEPT_MAP.get(c, "other") for c in combo]
        n_distinct = len(set(concepts))
        # Mean discrimination_hist (computed inline from hist_df)
        discriminations = []
        for c in combo:
            v = hist_df[c].to_numpy(dtype=float)
            q25, q75 = np.percentile(v, [25, 75])
            discriminations.append((q75 - q25) / _DSTAR)
        mean_discrim = float(np.mean(discriminations))
        composite = float(np.sqrt(max(mean_discrim, 0.0) * max(1.0 - median_rho, 0.0)))
        rows.append({
            "metrics": list(combo),
            "concepts": concepts,
            "median_abs_rho": median_rho,
            "max_abs_rho": max_rho,
            "n_distinct_concepts": n_distinct,
            "mean_discrimination": mean_discrim,
            "composite_score": composite,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--in-dir",
        type=Path,
        default=stage_output_dir(
            "02_calibration", "recompute_bounded_candidates", create=False,
        ),
        help="Directory containing recompute_bounded_candidates outputs.",
    )
    parser.add_argument(
        "--K",
        type=int,
        default=4,
        help="K-set cardinality (default: 4).",
    )
    parser.add_argument(
        "--gate-discrim",
        type=float,
        default=0.20,
        help="Per-candidate gate: minimum discrimination_hist (IQR/D*).",
    )
    parser.add_argument(
        "--gate-tail-res",
        type=float,
        default=1e-4,
        help="Per-candidate gate: minimum tail_resolution.",
    )
    parser.add_argument(
        "--gate-flood",
        type=float,
        default=0.05,
        help="Per-candidate gate: maximum flood_corner_frac.",
    )
    parser.add_argument(
        "--kset-rho-cap",
        type=float,
        default=0.60,
        help="K-set filter: max acceptable median |Spearman ρ|.",
    )
    parser.add_argument(
        "--kset-max-rho-cap",
        type=float,
        default=0.55,
        help=(
            "K-set filter: max acceptable MAX |Spearman ρ| within the set. "
            "Tighter than the median cap because a single redundant pair "
            "renders one of K objectives effectively non-informative in "
            "rank space. Default 0.55."
        ),
    )
    parser.add_argument(
        "--kset-min-concepts",
        type=int,
        default=4,
        help="K-set filter: minimum number of distinct concept tags.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Persist top-N K-sets per mapping in CSV outputs.",
    )
    args = parser.parse_args()

    out_dir = stage_output_dir("02_calibration", "select_bounded_kset")
    print(f"[select_bounded_kset] reading from {args.in_dir}")
    hist_df, archive_dfs = _load_recompute_outputs(args.in_dir)
    print(f"  {len(hist_df)} historical years, "
          f"{sum(len(d) for d in archive_dfs.values())} total Pareto rows "
          f"across {len(archive_dfs)} archives")

    print("[select_bounded_kset] computing per-candidate scores ...")
    cand_df = _per_candidate_scores(hist_df, archive_dfs)
    cand_path = out_dir / "per_candidate_scores.csv"
    cand_df.to_csv(cand_path, index=False)
    print(f"  wrote {cand_path}")

    # Apply per-candidate gates
    survivors = cand_df[
        (cand_df["discrimination_hist"] >= args.gate_discrim)
        & (cand_df["tail_resolution"] >= args.gate_tail_res)
        & (cand_df["flood_corner_frac"] <= args.gate_flood)
    ]
    n_survivors = len(survivors)
    print(f"[select_bounded_kset] {n_survivors}/{len(cand_df)} candidates "
          f"survive per-candidate gates")
    if n_survivors == 0:
        print("[select_bounded_kset] WARNING: no survivors. "
              "Loosen gates or expand the candidate pool. "
              "Writing empty selection_report.")
        (out_dir / "selection_report.json").write_text(json.dumps({
            "winner": None,
            "n_survivors": 0,
            "gates": vars(args),
        }, indent=2, default=str))
        return

    overall_results: Dict[str, Dict] = {}
    for mapping in ("g", "e"):
        sub = survivors[survivors["mapping"] == mapping]
        cand_list = list(sub["candidate"])
        print(f"[select_bounded_kset] enumerating K=4 subsets for mapping "
              f"'{mapping}' over {len(cand_list)} survivors ...")
        ksets = _enumerate_ksets(cand_list, hist_df, K=args.K)
        if ksets.empty:
            print(f"  no K=4 subsets producible for mapping '{mapping}'")
            overall_results[mapping] = {"winner": None, "n_ksets": 0}
            continue
        # Filter and rank
        keep = (
            (ksets["median_abs_rho"] <= args.kset_rho_cap)
            & (ksets["max_abs_rho"] <= args.kset_max_rho_cap)
            & (ksets["n_distinct_concepts"] >= args.kset_min_concepts)
        )
        passed = ksets[keep].sort_values("composite_score", ascending=False)
        print(f"  enumerated {len(ksets)} subsets; "
              f"{len(passed)} pass median_rho ≤ {args.kset_rho_cap}, "
              f"max_rho ≤ {args.kset_max_rho_cap}, "
              f"concepts ≥ {args.kset_min_concepts}")
        ranking_path = out_dir / f"kset_rankings_{mapping}.csv"
        # Convert list-valued columns to comma-joined strings for CSV legibility
        export = passed.head(args.top_n).copy()
        export["metrics"] = export["metrics"].apply(lambda xs: ",".join(xs))
        export["concepts"] = export["concepts"].apply(lambda xs: ",".join(xs))
        export.to_csv(ranking_path, index=False)
        print(f"  wrote {ranking_path}")

        if not passed.empty:
            best = passed.iloc[0]
            overall_results[mapping] = {
                "n_ksets_total": int(len(ksets)),
                "n_ksets_passed": int(len(passed)),
                "winner": {
                    "metrics": list(best["metrics"]),
                    "concepts": list(best["concepts"]),
                    "median_abs_rho": float(best["median_abs_rho"]),
                    "max_abs_rho": float(best["max_abs_rho"]),
                    "mean_discrimination": float(best["mean_discrimination"]),
                    "composite_score": float(best["composite_score"]),
                },
            }
        else:
            overall_results[mapping] = {
                "n_ksets_total": int(len(ksets)),
                "n_ksets_passed": 0,
                "winner": None,
            }

    # Compare mappings: pick the one with higher composite_score; tie-break on Mapping G (simpler).
    g_score = overall_results.get("g", {}).get("winner", {}).get("composite_score", -1.0) \
        if overall_results.get("g", {}).get("winner") else -1.0
    e_score = overall_results.get("e", {}).get("winner", {}).get("composite_score", -1.0) \
        if overall_results.get("e", {}).get("winner") else -1.0
    if g_score < 0 and e_score < 0:
        overall_winner = None
    elif g_score >= e_score * 0.95:
        overall_winner = {"mapping": "g", "details": overall_results["g"]["winner"]}
    else:
        overall_winner = {"mapping": "e", "details": overall_results["e"]["winner"]}

    report = {
        "K": args.K,
        "gates": {
            "discrim": args.gate_discrim,
            "tail_res": args.gate_tail_res,
            "flood": args.gate_flood,
            "kset_rho_cap": args.kset_rho_cap,
            "kset_max_rho_cap": args.kset_max_rho_cap,
            "kset_min_concepts": args.kset_min_concepts,
        },
        "n_candidates": len(cand_df),
        "n_survivors": int(n_survivors),
        "per_mapping": overall_results,
        "overall_winner": overall_winner,
    }
    report_path = out_dir / "selection_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str))
    print(f"[select_bounded_kset] wrote {report_path}")

    if overall_winner is not None:
        d = overall_winner["details"]
        print()
        print(f"=== OVERALL WINNER: mapping {overall_winner['mapping']!r} ===")
        print(f"  metrics: {d['metrics']}")
        print(f"  concepts: {d['concepts']}")
        print(f"  composite_score = {d['composite_score']:.4f}")
        print(f"  median |ρ| = {d['median_abs_rho']:.3f}; "
              f"max |ρ| = {d['max_abs_rho']:.3f}")
        print(f"  mean discrimination = {d['mean_discrimination']:.3f}")
    else:
        print("[select_bounded_kset] NO WINNER found — escalate to user.")


if __name__ == "__main__":
    main()
