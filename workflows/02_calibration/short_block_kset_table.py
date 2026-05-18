"""DD-15b — short-block K-set candidate table at T ∈ {1, 2}.

Enumerates K=3 and K=4 candidate metric tuples on the Kirsch correlation
matrix at each T, ranks them by the same composite score used in the
T=5 decision matrix (min spread × concept-diversity × low-correlation),
and writes a summary CSV. No SLURM compute beyond a few seconds.

Output:
  outputs/02_calibration/short_block_kset_table/short_block_kset_table.csv
  outputs/02_calibration/short_block_kset_table/short_block_kset_recommendation.json
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_CONCEPT_MAP,
    SHORT_BLOCK_METRIC_NAMES,
)


def _hist(T: int) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "02_calibration", "short_block_screening", slug=f"T{T}",
        create=False,
    ) / "block_chars_short.csv"
    return pd.read_csv(p) if p.exists() else None


def _kirsch(T: int, n_traces: int = 10_000, seed: int = 42) -> Optional[pd.DataFrame]:
    p = stage_output_dir(
        "03_kirsch_library", "build_short_block_library",
        slug=f"n{n_traces}_s{seed}", create=False,
    ) / f"characteristics_short_T{T}.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return pd.DataFrame(z["values"], columns=[str(n) for n in z["metric_names"]])


def _spread_score(values: np.ndarray) -> float:
    if values.size < 2:
        return 0.0
    med = float(np.median(values))
    std = float(np.std(values, ddof=1))
    q1, q3 = np.percentile(values, [25.0, 75.0])
    iqr = float(q3 - q1)
    denom = abs(med) + std
    return float(iqr / denom) if denom > 1e-12 else 0.0


def _enumerate(rho: pd.DataFrame, hist_chars: pd.DataFrame, K: int,
               cap: float, exclude_metrics: Optional[set] = None) -> List[Dict]:
    metrics = [m for m in rho.columns
               if (exclude_metrics is None or m not in exclude_metrics)]
    sub_rho = rho.loc[metrics, metrics]
    spread_lookup = {
        m: _spread_score(hist_chars[m].astype(float).dropna().values)
        for m in metrics if m in hist_chars.columns
    }
    concept_lookup = {m: SHORT_BLOCK_CONCEPT_MAP.get(m, "?") for m in metrics}

    cands: List[Dict] = []
    abs_rho = sub_rho.abs().values.copy()
    np.fill_diagonal(abs_rho, 0.0)
    name_to_idx = {m: i for i, m in enumerate(metrics)}
    for combo in combinations(metrics, K):
        # All distinct concepts (strict rung)
        if len({concept_lookup[m] for m in combo}) < K:
            continue
        idxs = [name_to_idx[m] for m in combo]
        max_off = 0.0
        for i in range(K):
            for j in range(i + 1, K):
                v = abs_rho[idxs[i], idxs[j]]
                if v > max_off:
                    max_off = v
        if max_off >= cap:
            continue
        spreads = [spread_lookup.get(m, 0.0) for m in combo]
        cands.append({
            "metrics": list(combo),
            "concepts": [concept_lookup[m] for m in combo],
            "spreads": spreads,
            "min_spread": float(min(spreads)),
            "sum_spread": float(sum(spreads)),
            "max_pairwise_rho": float(max_off),
        })
    cands.sort(key=lambda c: (-c["min_spread"], -c["sum_spread"]))
    return cands


def _composite(c: Dict, max_min_spread: float) -> float:
    score_spread = max(c["min_spread"] / max_min_spread, 1e-3) \
        if max_min_spread > 0 else 0.0
    score_indep = max(1.0 - c["max_pairwise_rho"], 1e-3)
    return float(np.exp(np.log([score_spread, score_indep]).mean()))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cap", type=float, nargs="+",
                   default=[0.4, 0.5, 0.6, 0.7])
    p.add_argument("--T-grid", type=int, nargs="+", default=[1, 2])
    p.add_argument("--top-N", type=int, default=10)
    p.add_argument("--exclude-T1-redundant", action="store_true",
                   default=True,
                   help="Drop total_flow_neg at T=1 (perfectly correlated "
                        "with min_annual_zscore there).")
    args = p.parse_args()

    out_dir = stage_output_dir(
        "02_calibration", "short_block_kset_table",
    )

    rows: List[Dict] = []
    rec_per_setting: Dict = {}
    for T in args.T_grid:
        hist = _hist(T)
        kir = _kirsch(T)
        if hist is None or kir is None:
            print(f"[warn] T={T}: missing artifacts; skipping")
            continue
        cols = [m for m in SHORT_BLOCK_METRIC_NAMES
                if m in kir.columns and m in hist.columns]
        rho = kir[cols].astype(float).corr(method="spearman")

        exclude = set()
        if args.exclude_T1_redundant and T == 1:
            exclude.add("total_flow_neg")  # perfectly redundant with min_annual_zscore

        for cap in args.cap:
            for K in (3, 4):
                cands = _enumerate(rho, hist, K=K, cap=cap,
                                   exclude_metrics=exclude)
                if not cands:
                    rows.append({
                        "T_years": T, "K": K, "cap": cap, "rank": None,
                        "metrics": "", "concepts": "",
                        "min_spread": np.nan,
                        "max_pairwise_rho": np.nan,
                        "composite_score": np.nan,
                    })
                    continue
                max_ms = max(c["min_spread"] for c in cands)
                for r, c in enumerate(cands[: args.top_N], 1):
                    rows.append({
                        "T_years": T, "K": K, "cap": cap, "rank": r,
                        "metrics": "|".join(c["metrics"]),
                        "concepts": "|".join(c["concepts"]),
                        "min_spread": c["min_spread"],
                        "sum_spread": c["sum_spread"],
                        "max_pairwise_rho": c["max_pairwise_rho"],
                        "composite_score": _composite(c, max_ms),
                    })

        # Per-T best at the highest-cap that admits both K=3 and K=4:
        # find the smallest cap that admits both (then pick best by score
        # at that cap).
        best_T: Dict = {"T_years": T, "k3": None, "k4": None}
        for cap in sorted(args.cap):
            for K in (3, 4):
                cands = _enumerate(rho, hist, K=K, cap=cap,
                                   exclude_metrics=exclude)
                if cands and best_T[f"k{K}"] is None:
                    max_ms = max(c["min_spread"] for c in cands)
                    top = max(cands, key=lambda c: _composite(c, max_ms))
                    best_T[f"k{K}"] = {
                        **top, "cap": cap,
                        "composite_score": _composite(top, max_ms),
                    }
        rec_per_setting[f"T{T}"] = best_T

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "short_block_kset_table.csv", index=False)
    (out_dir / "short_block_kset_recommendation.json").write_text(
        json.dumps(rec_per_setting, indent=2, default=lambda o: float(o)
                   if isinstance(o, np.floating) else None)
    )
    print(f"[short_kset] wrote {out_dir / 'short_block_kset_table.csv'} "
          f"({len(df)} rows)")
    for setting, rec in rec_per_setting.items():
        print(f"  {setting}:")
        for K_label, top in rec.items():
            if top is None or K_label == "T_years":
                continue
            if isinstance(top, dict):
                print(f"    {K_label} (cap={top.get('cap')}, "
                      f"composite={top.get('composite_score'):.3f}): "
                      f"{top.get('metrics')}")


if __name__ == "__main__":
    main()
