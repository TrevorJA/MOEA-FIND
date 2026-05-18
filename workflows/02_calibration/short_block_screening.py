"""DD-15b — historical short-block screening at T ∈ {1, 2}.

Companion to ``t_sensitivity_historical.py`` for the user-driven
sub-annual / single-year reframe (2026-04-29). Differences from the
T=5 pipeline:

* Metric pool is :data:`src.metrics.short_block.SHORT_BLOCK_METRIC_NAMES`
  (raw-flow Tier H + SSI-3 single-block Tier I), not the 28-metric
  multi-event candidate library.
* Each historical block is extracted with **3 burn-in months** prepended
  so the SSI-3 calculator (fitted on the full record) has a valid
  3-month rolling reference at every evaluation month.
* T=1 stride-1 yields 72 disjoint water years (each year used exactly
  once across blocks) — no autocorrelation issue.
* T=2 stride-1 yields 71 overlapping blocks; T=2 stride-2 yields 36
  disjoint blocks per tiling × 2 tilings = 72 unique observations.

Outputs (under ``outputs/02_calibration/short_block_screening/T{T}/``):

* ``block_chars_short.csv`` — per-block metric matrix
* ``per_metric_spread.csv`` — robust spread + screen flag
* ``spearman_corr.csv``, ``pearson_corr.csv`` — overlapping-block
  correlation matrices
* ``spearman_corr_disjoint.csv`` — disjoint-tiling-mean correlations
  (statistically unbiased reference)
* ``config.json`` — reproducibility metadata
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.extended import FullRecordRefs  # noqa: E402
from src.metrics.screening import per_metric_spread  # noqa: E402
from src.metrics.objectives import flows_to_series, make_ssi_calculator  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.metrics.short_block import (  # noqa: E402
    SHORT_BLOCK_CONCEPT_MAP,
    SHORT_BLOCK_METRIC_NAMES,
    compute_short_block_metrics,
)

STAGE = "02_calibration"
DRIVER = "short_block_screening"
BURNIN_MONTHS = 3


def _historical_overlap_blocks(
    monthly_1d: np.ndarray, T_years: int, burnin_months: int = BURNIN_MONTHS,
):
    """Yield stride-1 (block_eval, ssi_window) tuples, skipping blocks
    where the 3-month burn-in would precede the start of the record."""
    n_years = len(monthly_1d) // 12
    for j in range(1, n_years - T_years + 1):
        eval_start = j * 12
        eval_end = eval_start + T_years * 12
        burnin_start = eval_start - burnin_months
        if burnin_start < 0:
            continue
        yield (
            j,
            monthly_1d[eval_start:eval_end].copy(),
            monthly_1d[burnin_start:eval_end].copy(),
        )


def _historical_disjoint_blocks(
    monthly_1d: np.ndarray, T_years: int, burnin_months: int = BURNIN_MONTHS,
):
    """Yield disjoint stride-T tilings (each year used at most once per
    tiling) as a list of lists of (block_eval, ssi_window) tuples."""
    n_years = len(monthly_1d) // 12
    tilings: List[List[tuple]] = []
    for offset in range(T_years):
        tiling: List[tuple] = []
        j = max(1, offset)  # need ≥1 year of preceding data for burn-in
        if j == 0:
            j = T_years  # skip first tiling pos for offset=0 too
        # Step in T-year increments
        while j + T_years <= n_years:
            eval_start = j * 12
            eval_end = eval_start + T_years * 12
            burnin_start = eval_start - burnin_months
            if burnin_start < 0:
                j += T_years
                continue
            tiling.append((
                j,
                monthly_1d[eval_start:eval_end].copy(),
                monthly_1d[burnin_start:eval_end].copy(),
            ))
            j += T_years
        if tiling:
            tilings.append(tiling)
    return tilings


def _compute_block_matrix(blocks_with_burnin, ssi3_calc, refs) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    for j, block_eval, ssi_input in blocks_with_burnin:
        ssi_series = ssi3_calc.transform(
            flows_to_series(ssi_input, start_date="2100-01-01")
        )
        chars = compute_short_block_metrics(
            block_eval, ssi_series, refs,
            eval_first_idx_in_ssi=BURNIN_MONTHS,
        )
        chars["water_year_index"] = int(j)
        rows.append(chars)
    return pd.DataFrame(rows)


def _spearman_avg_over_tilings(tilings_block_chars: List[pd.DataFrame],
                                cols: List[str]) -> pd.DataFrame:
    if not tilings_block_chars:
        return pd.DataFrame()
    rho_sum = None
    n = 0
    for chars in tilings_block_chars:
        sub = chars[cols].astype(float)
        if len(sub) < 3:
            continue
        rho = sub.corr(method="spearman")
        rho_sum = rho if rho_sum is None else (rho_sum + rho)
        n += 1
    if n == 0:
        return pd.DataFrame()
    return rho_sum / n


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, required=True, choices=[1, 2])
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER, slug=f"T{args.T_years}")
    print(f"[02/{DRIVER}] T={args.T_years} → {out_dir}")

    cache = PROJECT_ROOT / "outputs" / "data_cache"
    cache.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache)

    ssi3 = make_ssi_calculator(timescale=3)
    ssi3.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    refs = FullRecordRefs.from_full_record(monthly_1d)

    # --- Stride-1 overlapping blocks (max sample) ---
    overlap_blocks = list(_historical_overlap_blocks(monthly_1d, args.T_years))
    chars_overlap = _compute_block_matrix(overlap_blocks, ssi3, refs)
    chars_overlap.to_csv(out_dir / "block_chars_short.csv", index=False)
    print(f"[diag] overlap: {len(chars_overlap)} blocks × "
          f"{len(SHORT_BLOCK_METRIC_NAMES)} metrics")

    # Per-metric spread + screen on the stride-1 sample.
    cols = [m for m in SHORT_BLOCK_METRIC_NAMES if m in chars_overlap.columns]
    # Inject the short-block CONCEPT_MAP into per_metric_spread by
    # patching the lookup; per_metric_spread imports CONCEPT_MAP from
    # extended_drought_metrics so we replicate the spread logic locally
    # for the new metric pool.
    spread_rows: List[Dict] = []
    from scipy import stats as scstats
    MIN_STD = 1e-9
    SKEW_DROP_THRESHOLD = 3.0
    for name in cols:
        v = chars_overlap[name].astype(float).values
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        mean = float(np.mean(v))
        med = float(np.median(v))
        std = float(np.std(v, ddof=1)) if v.size > 1 else 0.0
        q1, q3 = np.percentile(v, [25.0, 75.0])
        iqr = float(q3 - q1)
        denom = abs(med) + std
        spread_score = float(iqr / denom) if denom > 1e-12 else 0.0
        skew = float(scstats.skew(v, bias=False)) if v.size > 2 else 0.0
        passes = (
            std >= MIN_STD
            and iqr > 0.0
            and spread_score > 0.0
            and np.isfinite(spread_score)
            and abs(skew) <= SKEW_DROP_THRESHOLD
        )
        spread_rows.append({
            "metric": name,
            "concept": SHORT_BLOCK_CONCEPT_MAP.get(name, "unknown"),
            "mean": mean, "median": med, "std": std,
            "min": float(np.min(v)), "max": float(np.max(v)),
            "range": float(np.max(v) - np.min(v)),
            "iqr": iqr, "spread_score": spread_score, "skew": skew,
            "passes_screen": bool(passes),
        })
    spread_df = pd.DataFrame(spread_rows).sort_values(
        "spread_score", ascending=False, na_position="last"
    ).reset_index(drop=True)
    spread_df.to_csv(out_dir / "per_metric_spread.csv", index=False)
    print(f"[diag] {spread_df.passes_screen.sum()}/{len(spread_df)} "
          f"metrics pass screen")

    # --- Correlations on overlapping blocks ---
    sub = chars_overlap[cols].astype(float)
    pearson = sub.corr(method="pearson")
    spearman = sub.corr(method="spearman")
    pearson.to_csv(out_dir / "pearson_corr.csv")
    spearman.to_csv(out_dir / "spearman_corr.csv")

    # --- Disjoint-tiling correlations (statistically unbiased) ---
    if args.T_years > 1:
        tilings = _historical_disjoint_blocks(monthly_1d, args.T_years)
        per_tiling_chars = [_compute_block_matrix(t, ssi3, refs) for t in tilings]
        spearman_disj = _spearman_avg_over_tilings(per_tiling_chars, cols)
        if not spearman_disj.empty:
            spearman_disj.to_csv(out_dir / "spearman_corr_disjoint.csv")
            print(f"[diag] disjoint correlations: {len(tilings)} tilings × "
                  f"{len(per_tiling_chars[0]) if per_tiling_chars else 0} blocks")
    else:
        # T=1: stride-1 already gives disjoint independent blocks (each
        # water year used exactly once), so spearman == spearman_disjoint.
        spearman.to_csv(out_dir / "spearman_corr_disjoint.csv")
        print("[diag] T=1: stride-1 blocks are already disjoint")

    # --- Config ---
    cfg = {
        "T_years": int(args.T_years),
        "burnin_months": BURNIN_MONTHS,
        "n_blocks_overlap": int(len(chars_overlap)),
        "metric_pool_size": len(cols),
        "metric_names": cols,
        "concept_map": {m: SHORT_BLOCK_CONCEPT_MAP.get(m, "unknown") for m in cols},
        "n_pass_screen": int(spread_df.passes_screen.sum()),
    }
    (out_dir / "config.json").write_text(json.dumps(cfg, indent=2))
    print(f"[diag] done. wrote {out_dir / 'config.json'}")


if __name__ == "__main__":
    main()
