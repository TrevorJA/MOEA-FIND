"""Stage 06 / verify_drought_coverage -- pre-flight Pareto verification.

Implements Criteria 1-5 from docs/moea_find_verification_criteria.md against
a completed MOEA-FIND archive (results.json produced by stage 04). Produces a
single ``verification_report.json`` plus the auxiliary numeric arrays needed
by the paired plotting driver under
``src/plotting/06_pywrdrb_reeval/verify_drought_coverage.py``.

This driver writes only numerical artifacts -- no matplotlib calls.

Subsets:
    drought_subset = union of (D_j >= historical median) across K objectives
    nominal_subset = strict intersection complement (all D_j < historical median)

Each criterion reports ``pass`` / ``fail`` / ``insufficient_data`` plus the
underlying numeric diagnostics so the report can be diffed across runs.

Outputs under ``outputs/06_pywrdrb_reeval/verify_drought_coverage/<src_slug>/``:
    - config.json
    - verification_report.json
    - subsets.npz             (drought_mask, nominal_mask, hist_medians)
    - annual_means.npz        (A_syn, A_hist for the spread plot)
    - hist_block_chars.npz    (cached for downstream reuse)

Usage:
    python workflows/06_pywrdrb_reeval/verify_drought_coverage.py \\
        --pareto-results outputs/04_moea_find_single_site/run_moea_find/<slug>/results.json
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.hydrology.historical_blocks import (  # noqa: E402
    resample_historical_blocks,
    resample_historical_blocks_2d,
    compute_historical_block_chars,
)
from src.experiment import prepare_data  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    compute_ssi,
    flows_to_series,
    get_drought_metrics,
)
from src.optimization.constraints import _annual_totals, _lag1_ac  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "06_pywrdrb_reeval"
DRIVER = "verify_drought_coverage"

EXCEEDANCE_LOW = (0.70, 0.80, 0.90, 0.95, 0.99)
EXCEEDANCE_HIGH = (0.01, 0.05, 0.10, 0.20, 0.30)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _fdc_at_exceedance(trace: np.ndarray, exceedances: Sequence[float]) -> np.ndarray:
    """Return flow values at the requested exceedance probabilities.

    Uses the empirical grid ``i / (n+1)`` for ``i = 1..n`` and linear
    interpolation; this matches :func:`src.plotting.trace_diagnostics._fdc`.
    """
    sorted_desc = np.sort(trace)[::-1]
    n = len(sorted_desc)
    grid = np.arange(1, n + 1) / (n + 1)
    return np.interp(np.asarray(exceedances, dtype=float), grid, sorted_desc)


def _percentile_safe(arr: np.ndarray, q: float) -> float:
    if arr.size == 0:
        return float("nan")
    return float(np.percentile(arr, q))


def _bbox_overlap_ratio(syn: np.ndarray, hist: np.ndarray) -> float:
    """Volume of (syn_bbox intersect hist_bbox) / volume of hist_bbox."""
    lo = np.maximum(syn.min(axis=0), hist.min(axis=0))
    hi = np.minimum(syn.max(axis=0), hist.max(axis=0))
    inter = np.clip(hi - lo, 0.0, None)
    hist_span = np.clip(hist.max(axis=0) - hist.min(axis=0), 1e-12, None)
    ratio = float(np.prod(inter / hist_span))
    return ratio


def _ssi_onset_months(
    trace_1d: np.ndarray,
    ssi_calc,
    start_date: str,
    end_drought_threshold_months: int = 3,
) -> np.ndarray:
    """Return a 1D array of drought-event start calendar-months (1-12)."""
    ssi = ssi_calc.transform(flows_to_series(trace_1d, start_date=start_date))
    dm = get_drought_metrics(
        ssi, end_drought_threshold_months=end_drought_threshold_months,
    )
    if len(dm) == 0 or "start" not in dm.columns:
        return np.empty(0, dtype=int)
    return pd.to_datetime(dm["start"]).dt.month.values.astype(int)


# ---------------------------------------------------------------------------
# Criterion implementations
# ---------------------------------------------------------------------------


def criterion_1_coverage(
    pareto_drought_metrics: np.ndarray,
    hist_block_chars: np.ndarray,
    objective_keys: Sequence[str],
) -> Dict:
    """Drought-space coverage (Criterion 1)."""
    ratio = _bbox_overlap_ratio(pareto_drought_metrics, hist_block_chars)
    p95 = np.percentile(hist_block_chars, 95, axis=0)
    beyond_p95 = (pareto_drought_metrics >= p95).any(axis=0)
    syn_range = pareto_drought_metrics.max(axis=0) - pareto_drought_metrics.min(axis=0)
    hist_range = hist_block_chars.max(axis=0) - hist_block_chars.min(axis=0)
    range_ratio = syn_range / np.clip(hist_range, 1e-12, None)

    per_objective = []
    for j, key in enumerate(objective_keys):
        per_objective.append({
            "objective": key,
            "hist_p95": float(p95[j]),
            "syn_max": float(pareto_drought_metrics[:, j].max()),
            "beyond_hist_p95": bool(beyond_p95[j]),
            "syn_range": float(syn_range[j]),
            "hist_range": float(hist_range[j]),
            "range_ratio": float(range_ratio[j]),
            "range_ratio_ge_2": bool(range_ratio[j] >= 2.0),
        })

    passes = (
        ratio > 0.50
        and all(bool(b) for b in beyond_p95)
        and all(bool(r >= 2.0) for r in range_ratio)
    )
    return {
        "criterion": 1,
        "name": "drought_space_coverage",
        "status": "pass" if passes else "fail",
        "bbox_overlap_ratio": float(ratio),
        "bbox_overlap_threshold": 0.50,
        "per_objective": per_objective,
    }


def criterion_2_low_flow_directionality(
    drought_traces_1d: List[np.ndarray],
    hist_blocks_1d: List[np.ndarray],
) -> Dict:
    """Drought subset low-flow directionality (Criterion 2)."""
    if len(drought_traces_1d) < 100:
        return {
            "criterion": 2,
            "name": "low_flow_directionality",
            "status": "insufficient_data",
            "n_drought_subset": len(drought_traces_1d),
            "message": "fewer than 100 drought_subset traces; ramp NFE",
        }

    syn_fdc = np.array([
        _fdc_at_exceedance(t, EXCEEDANCE_LOW) for t in drought_traces_1d
    ])
    hist_fdc = np.array([
        _fdc_at_exceedance(b, EXCEEDANCE_LOW) for b in hist_blocks_1d
    ])

    syn_p50 = np.percentile(syn_fdc, 50, axis=0)
    syn_p10 = np.percentile(syn_fdc, 10, axis=0)
    hist_p50 = np.percentile(hist_fdc, 50, axis=0)
    hist_p10 = np.percentile(hist_fdc, 10, axis=0)

    per_exceedance = []
    all_pass = True
    for j, e in enumerate(EXCEEDANCE_LOW):
        ok = bool(syn_p50[j] <= hist_p50[j] and syn_p10[j] <= hist_p10[j])
        all_pass = all_pass and ok
        per_exceedance.append({
            "exceedance": float(e),
            "syn_p50": float(syn_p50[j]),
            "hist_p50": float(hist_p50[j]),
            "syn_p10": float(syn_p10[j]),
            "hist_p10": float(hist_p10[j]),
            "pass": ok,
        })

    return {
        "criterion": 2,
        "name": "low_flow_directionality",
        "status": "pass" if all_pass else "fail",
        "n_drought_subset": len(drought_traces_1d),
        "per_exceedance": per_exceedance,
    }


def criterion_3_nominal_fidelity(
    nominal_traces_1d: List[np.ndarray],
    nominal_traces_2d: List[np.ndarray],
    hist_blocks_1d: List[np.ndarray],
    hist_blocks_2d: List[np.ndarray],
) -> Dict:
    """Nominal-subset non-drought hydrologic fidelity (Criterion 3)."""
    if len(nominal_traces_1d) < 100:
        return {
            "criterion": 3,
            "name": "non_drought_fidelity",
            "status": "insufficient_data",
            "n_nominal_subset": len(nominal_traces_1d),
            "message": "fewer than 100 nominal_subset traces; ramp NFE",
        }

    # ---- High-flow FDC envelope containment (>= 90% of traces inside band) ----
    syn_fdc = np.array([
        _fdc_at_exceedance(t, EXCEEDANCE_HIGH) for t in nominal_traces_1d
    ])
    hist_fdc = np.array([
        _fdc_at_exceedance(b, EXCEEDANCE_HIGH) for b in hist_blocks_1d
    ])
    hist_lo = np.percentile(hist_fdc, 10, axis=0)
    hist_hi = np.percentile(hist_fdc, 90, axis=0)
    inside = (syn_fdc >= hist_lo) & (syn_fdc <= hist_hi)
    pct_inside_per_excd = inside.mean(axis=0) * 100.0

    fdc_pass = all(p >= 90.0 for p in pct_inside_per_excd)

    # ---- Lag-1 AC intervals overlap ----
    syn_ac = np.array([_lag1_ac(t) for t in nominal_traces_1d])
    hist_ac = np.array([_lag1_ac(b) for b in hist_blocks_1d])
    syn_lo, syn_hi = _percentile_safe(syn_ac, 5), _percentile_safe(syn_ac, 95)
    hist_ac_lo, hist_ac_hi = _percentile_safe(hist_ac, 5), _percentile_safe(hist_ac, 95)
    ac_pass = bool(syn_lo <= hist_ac_hi and syn_hi >= hist_ac_lo)

    # ---- Seasonal cycle containment (mean and std per calendar month) ----
    syn_month_mean = np.array([t.mean(axis=0) for t in nominal_traces_2d])  # (n, 12)
    syn_month_std = np.array([t.std(axis=0, ddof=1) for t in nominal_traces_2d])
    hist_month_mean = np.array([b.mean(axis=0) for b in hist_blocks_2d])
    hist_month_std = np.array([b.std(axis=0, ddof=1) for b in hist_blocks_2d])

    def envelope_containment(syn: np.ndarray, hist: np.ndarray) -> float:
        lo = np.percentile(hist, 10, axis=0)
        hi = np.percentile(hist, 90, axis=0)
        inside_ = (syn >= lo) & (syn <= hi)
        return float(inside_.mean(axis=0).min() * 100.0)  # worst month

    mean_inside_pct = envelope_containment(syn_month_mean, hist_month_mean)
    std_inside_pct = envelope_containment(syn_month_std, hist_month_std)
    seasonal_pass = mean_inside_pct >= 90.0 and std_inside_pct >= 90.0

    all_pass = bool(fdc_pass and ac_pass and seasonal_pass)
    return {
        "criterion": 3,
        "name": "non_drought_fidelity",
        "status": "pass" if all_pass else "fail",
        "n_nominal_subset": len(nominal_traces_1d),
        "fdc_high_flow": {
            "exceedances": list(EXCEEDANCE_HIGH),
            "pct_inside_per_exceedance": [float(p) for p in pct_inside_per_excd],
            "pass": fdc_pass,
        },
        "lag1_ac": {
            "syn_5_95": [syn_lo, syn_hi],
            "hist_5_95": [hist_ac_lo, hist_ac_hi],
            "intervals_overlap": ac_pass,
        },
        "seasonal_cycle": {
            "min_month_pct_inside_mean": float(mean_inside_pct),
            "min_month_pct_inside_std": float(std_inside_pct),
            "pass": seasonal_pass,
        },
    }


def criterion_4_flow_space_spread(
    pareto_traces_1d: List[np.ndarray],
    hist_blocks_1d: List[np.ndarray],
    drought_metrics: np.ndarray,
    hist_block_chars: np.ndarray,
    drought_mask: np.ndarray,
    nominal_mask: np.ndarray,
):
    """Flow-space spread + corner bifurcation (Criterion 4).

    Returns ``(report_dict, A_syn, A_hist)`` so the plotting driver can
    re-render the spread histogram from cached arrays without re-running
    ``_annual_totals``.
    """
    A_syn = np.array([_annual_totals(t).mean() for t in pareto_traces_1d])
    A_hist = np.array([_annual_totals(b).mean() for b in hist_blocks_1d])

    std_syn = float(A_syn.std(ddof=1)) if len(A_syn) > 1 else 0.0
    std_hist = float(A_hist.std(ddof=1)) if len(A_hist) > 1 else 0.0
    spread_ratio = std_syn / std_hist if std_hist > 0 else float("nan")
    spread_pass = bool(0.9 <= spread_ratio <= 1.1)

    hist_med_annual = float(np.median(A_hist))
    corner_diag = {
        "hist_median_annual_mean": hist_med_annual,
        "low_duration_corner_mean_annual": None,
        "high_duration_corner_mean_annual": None,
        "bifurcation_pass": None,
    }
    if drought_mask.any() and nominal_mask.any():
        low_corner_mean = float(A_syn[nominal_mask].mean())
        high_corner_mean = float(A_syn[drought_mask].mean())
        corner_diag["low_duration_corner_mean_annual"] = low_corner_mean
        corner_diag["high_duration_corner_mean_annual"] = high_corner_mean
        corner_diag["bifurcation_pass"] = bool(
            low_corner_mean >= hist_med_annual
            and high_corner_mean <= hist_med_annual
        )

    overall_pass = spread_pass and (
        corner_diag["bifurcation_pass"] in (True, None)
    )
    report = {
        "criterion": 4,
        "name": "flow_space_spread",
        "status": "pass" if overall_pass else "fail",
        "std_syn": std_syn,
        "std_hist": std_hist,
        "spread_ratio": spread_ratio,
        "spread_pass": spread_pass,
        "spread_tolerance": [0.9, 1.1],
        "corner_bifurcation": corner_diag,
    }
    return report, A_syn, A_hist


def criterion_5_drought_onset(
    drought_traces_1d: List[np.ndarray],
    hist_blocks_1d: List[np.ndarray],
    ssi_calc,
    start_date: str,
    end_drought_threshold_months: int = 3,
) -> Dict:
    """Drought-onset month chi-squared test (Criterion 5)."""
    from scipy.stats import chi2_contingency

    if len(drought_traces_1d) < 100:
        return {
            "criterion": 5,
            "name": "drought_onset_seasonality",
            "status": "insufficient_data",
            "n_drought_subset": len(drought_traces_1d),
        }

    syn_onsets = np.concatenate([
        _ssi_onset_months(t, ssi_calc, start_date,
                          end_drought_threshold_months)
        for t in drought_traces_1d
    ]) if drought_traces_1d else np.empty(0, dtype=int)
    hist_onsets = np.concatenate([
        _ssi_onset_months(b, ssi_calc, start_date,
                          end_drought_threshold_months)
        for b in hist_blocks_1d
    ]) if hist_blocks_1d else np.empty(0, dtype=int)

    if len(syn_onsets) < 12 or len(hist_onsets) < 12:
        return {
            "criterion": 5,
            "name": "drought_onset_seasonality",
            "status": "insufficient_data",
            "n_synthetic_events": int(len(syn_onsets)),
            "n_historical_events": int(len(hist_onsets)),
            "message": "too few drought events for chi-squared test",
        }

    syn_hist = np.bincount(syn_onsets, minlength=13)[1:13]  # months 1..12
    hist_hist = np.bincount(hist_onsets, minlength=13)[1:13]

    # Guard: chi2 requires no entirely-empty column. Add a tiny constant
    # then round to keep pytest semantics clean when real counts are zero.
    eps = 1e-9
    table = np.vstack([syn_hist, hist_hist]).astype(float) + eps
    chi2, p_value, dof, _ = chi2_contingency(table)

    return {
        "criterion": 5,
        "name": "drought_onset_seasonality",
        "status": "pass" if p_value >= 0.05 else "fail",
        "p_value": float(p_value),
        "chi2_statistic": float(chi2),
        "dof": int(dof),
        "n_synthetic_events": int(len(syn_onsets)),
        "n_historical_events": int(len(hist_onsets)),
        "synthetic_histogram": [int(x) for x in syn_hist],
        "historical_histogram": [int(x) for x in hist_hist],
    }


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pareto-results", type=Path, required=True,
                    help="Path to results.json from a MOEA-FIND stage 04 run.")
    ap.add_argument("--min-subset-size", type=int, default=100,
                    help="Minimum traces per subset for Criteria 2, 3, 5 "
                         "(default: 100).")
    args = ap.parse_args()

    if not args.pareto_results.exists():
        raise FileNotFoundError(args.pareto_results)
    src_slug = args.pareto_results.parent.name
    out_dir = stage_output_dir(STAGE, DRIVER, src_slug)

    print(f"[06/verify] loading {args.pareto_results}")
    results = json.loads(args.pareto_results.read_text())
    pareto_traces_1d = [np.asarray(t, dtype=float) for t in results["pareto_traces_1d"]]
    pareto_traces_2d = [np.asarray(t, dtype=float) for t in results["pareto_traces_2d"]]
    drought_metrics = np.asarray(results["drought_metrics"], dtype=float)
    objective_keys = tuple(results.get(
        "objective_keys", ["mean_duration", "mean_avg_severity"]))
    anti_ideal = np.asarray(results.get("anti_ideal", []), dtype=float)
    n_years = int(results["n_years_out"])
    ssi_timescale = int(results.get("ssi_timescale", 3))

    print(f"[06/verify] {len(pareto_traces_1d)} Pareto traces, "
          f"K={drought_metrics.shape[1]}, T={n_years} years, SSI-{ssi_timescale}")

    # Config dump
    config = {
        "stage": STAGE,
        "driver": DRIVER,
        "pareto_source": str(args.pareto_results),
        "src_slug": src_slug,
        "n_pareto": len(pareto_traces_1d),
        "n_years": n_years,
        "ssi_timescale": ssi_timescale,
        "objective_keys": list(objective_keys),
        "min_subset_size": args.min_subset_size,
    }
    (out_dir / "config.json").write_text(json.dumps(config, indent=2))

    # Historical data + blocks (match stage 04 conventions)
    print("[06/verify] loading historical flows + building T-year blocks")
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir=cache_dir)
    hist_blocks_1d = resample_historical_blocks(monthly_1d, T_years=n_years, stride=1)
    hist_blocks_2d = resample_historical_blocks_2d(monthly_2d, T_years=n_years, stride=1)
    print(f"[06/verify] {len(hist_blocks_1d)} historical blocks at T={n_years}")

    # Historical block drought chars -- load cache from the source archive if
    # present (stage 04 writes it next to results.json), else recompute.
    src_dir = args.pareto_results.parent
    hist_char_npz_src = src_dir / "historical_block_chars.npz"
    hist_block_chars = None
    ssi_calc = None
    if hist_char_npz_src.exists():
        print(f"[06/verify] loading cached {hist_char_npz_src}")
        npz = np.load(hist_char_npz_src, allow_pickle=True)
        hist_block_chars = np.asarray(npz["chars"], dtype=float)
        cached_keys = [str(k) for k in npz["objective_keys"]]
        if tuple(cached_keys) != objective_keys:
            warnings.warn(
                f"cached objective_keys {cached_keys} != results "
                f"{objective_keys}; recomputing"
            )
            hist_block_chars = None

    # Fit the same SSI calculator the MOEA used -- prefit on full history.
    _ssi_series, ssi_calc = compute_ssi(monthly_1d, timescale=ssi_timescale)

    if hist_block_chars is None:
        print("[06/verify] recomputing historical block drought chars")
        hist_block_chars = compute_historical_block_chars(
            monthly_1d, T_years=n_years, ssi_calc=ssi_calc,
            objective_keys=objective_keys, stride=1,
        )

    # Persist a local copy under our stage output for plot reuse.
    np.savez(
        out_dir / "hist_block_chars.npz",
        chars=hist_block_chars,
        objective_keys=np.array(list(objective_keys)),
    )

    # Subsets
    hist_medians = np.median(hist_block_chars, axis=0)
    drought_mask = (drought_metrics >= hist_medians).any(axis=1)
    nominal_mask = (drought_metrics < hist_medians).all(axis=1)
    n_drought = int(drought_mask.sum())
    n_nominal = int(nominal_mask.sum())
    print(f"[06/verify] subsets: drought={n_drought}, nominal={n_nominal} "
          f"(union vs intersection complement on hist medians={hist_medians})")

    drought_traces_1d = [pareto_traces_1d[i]
                         for i in np.where(drought_mask)[0]]
    nominal_traces_1d = [pareto_traces_1d[i]
                         for i in np.where(nominal_mask)[0]]
    nominal_traces_2d = [pareto_traces_2d[i]
                         for i in np.where(nominal_mask)[0]]

    # Criteria
    print("[06/verify] criterion 1 -- drought-space coverage")
    c1 = criterion_1_coverage(drought_metrics, hist_block_chars, objective_keys)

    print("[06/verify] criterion 2 -- low-flow directionality")
    c2 = criterion_2_low_flow_directionality(drought_traces_1d, hist_blocks_1d)

    print("[06/verify] criterion 3 -- non-drought fidelity")
    c3 = criterion_3_nominal_fidelity(
        nominal_traces_1d, nominal_traces_2d, hist_blocks_1d, hist_blocks_2d,
    )

    print("[06/verify] criterion 4 -- flow-space spread")
    c4, A_syn, A_hist = criterion_4_flow_space_spread(
        pareto_traces_1d, hist_blocks_1d, drought_metrics,
        hist_block_chars, drought_mask, nominal_mask,
    )

    print("[06/verify] criterion 5 -- drought onset seasonality (SSI re-transform)")
    c5 = criterion_5_drought_onset(
        drought_traces_1d, hist_blocks_1d, ssi_calc,
        start_date="2100-01-01",
    )

    # Persist arrays needed by the plotting driver.
    np.savez(
        out_dir / "subsets.npz",
        drought_mask=drought_mask,
        nominal_mask=nominal_mask,
        hist_medians=hist_medians,
        drought_metrics=drought_metrics,
        anti_ideal=anti_ideal,
    )
    np.savez(out_dir / "annual_means.npz", A_syn=A_syn, A_hist=A_hist)

    # Assemble report
    report = {
        "archive": str(args.pareto_results),
        "src_slug": src_slug,
        "n_pareto": len(pareto_traces_1d),
        "n_years": n_years,
        "ssi_timescale": ssi_timescale,
        "objective_keys": list(objective_keys),
        "historical_medians": [float(x) for x in hist_medians],
        "subsets": {
            "n_drought_subset": n_drought,
            "n_nominal_subset": n_nominal,
            "drought_rule": "union: any D_j >= historical median",
            "nominal_rule": "intersection complement: all D_j < historical median",
            "min_subset_size_target": args.min_subset_size,
        },
        "criteria": {
            "c1_coverage": c1,
            "c2_low_flow": c2,
            "c3_nominal_fidelity": c3,
            "c4_spread": c4,
            "c5_drought_onset": c5,
        },
    }
    overall_statuses = [c["status"] for c in report["criteria"].values()]
    report["overall_status"] = (
        "pass" if all(s == "pass" for s in overall_statuses)
        else ("fail" if "fail" in overall_statuses else "insufficient_data")
    )

    report_path = out_dir / "verification_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"[06/verify] wrote {report_path} -- overall={report['overall_status']}")

    # Compact terminal summary
    print("\n[06/verify] summary")
    for key, crit in report["criteria"].items():
        print(f"  {key:25s} {crit['status']}")
    print(f"  {'overall':25s} {report['overall_status']}")
    print(f"  outputs: {out_dir}")


if __name__ == "__main__":
    main()
