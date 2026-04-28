"""Stage A verification diagnostic — pass/fail report for MOEA-FIND archive.

Implements Criteria 1-5 from docs/moea_find_verification_criteria.md against
a completed MOEA-FIND archive (results.json produced by script 04 or 08).
Produces verification_report.json and five figures under the archive's
figures/ directory.

Subsets (per 2026-04-17 session notes in the verification criteria doc):
    drought_subset = union of (D_j >= historical median) across K objectives
    nominal_subset = strict intersection complement (all D_j < historical median)

Each criterion reports ``pass`` / ``fail`` / ``insufficient_data`` plus the
underlying numeric diagnostics so the report can be diffed across runs.

Usage:
    python workflows/06_pywrdrb_reeval/verify_drought_coverage.py \\
        outputs/exp04_kirsch_single_site/<variant>/results.json
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.historical_blocks import (  # noqa: E402
    resample_historical_blocks,
    resample_historical_blocks_2d,
    compute_historical_block_chars,
)
from src.experiment_utils import prepare_data  # noqa: E402
from src.objectives import (  # noqa: E402
    compute_ssi,
    flows_to_series,
    get_drought_metrics,
)
from src.constraints import _annual_totals, _lag1_ac  # noqa: E402

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
    """Volume of (syn_bbox ∩ hist_bbox) / volume of hist_bbox.

    ``syn`` / ``hist`` are ``(n_points, K)`` arrays. Non-cyclic only —
    cyclic axes (peak_severity_month) use wrapped overlap handled by
    the caller.
    """
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

    # ---- High-flow FDC envelope containment (≥90% of traces inside band) ----
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
) -> Dict:
    """Flow-space spread + corner bifurcation (Criterion 4)."""
    A_syn = np.array([_annual_totals(t).mean() for t in pareto_traces_1d])
    A_hist = np.array([_annual_totals(b).mean() for b in hist_blocks_1d])

    std_syn = float(A_syn.std(ddof=1)) if len(A_syn) > 1 else 0.0
    std_hist = float(A_hist.std(ddof=1)) if len(A_hist) > 1 else 0.0
    spread_ratio = std_syn / std_hist if std_hist > 0 else float("nan")
    spread_pass = bool(0.9 <= spread_ratio <= 1.1)

    # Corner bifurcation — use the *mean_duration* axis (first objective) as
    # the canonical drought-strength dimension. Nominal corner (wet, low D)
    # should be wetter than historical median; drought corner (high D) drier.
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
    return {
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
# Figures
# ---------------------------------------------------------------------------


def _save_fig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"[verify] wrote {path}")


def plot_c1_coverage(
    pareto_drought_metrics: np.ndarray,
    hist_block_chars: np.ndarray,
    objective_keys: Sequence[str],
    anti_ideal: np.ndarray,
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    K = pareto_drought_metrics.shape[1]
    n_panels = K * (K - 1) // 2 if K >= 2 else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5), squeeze=False)
    axes = axes.flatten()

    p = 0
    for i in range(K):
        for j in range(i + 1, K):
            ax = axes[p]
            ax.scatter(hist_block_chars[:, i], hist_block_chars[:, j],
                       s=40, c="#1f77b4", alpha=0.6, label="Historical blocks")
            ax.scatter(pareto_drought_metrics[:, i], pareto_drought_metrics[:, j],
                       s=8, c="#ff7f0e", alpha=0.5, label="Pareto archive")
            ax.scatter([anti_ideal[i]], [anti_ideal[j]],
                       marker="X", s=180, c="red", zorder=5, label="D*")
            ax.set_xlabel(objective_keys[i])
            ax.set_ylabel(objective_keys[j])
            ax.legend(fontsize=8, loc="best")
            p += 1
    fig.suptitle("Criterion 1 — Drought-space coverage", fontsize=12)
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c2_low_flow(
    report: Dict,
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    if report.get("status") == "insufficient_data":
        ax.text(0.5, 0.5, f"insufficient_data\nn_drought_subset={report.get('n_drought_subset', 0)}",
                ha="center", va="center", transform=ax.transAxes)
        _save_fig(fig, output_path)
        return

    xs = [r["exceedance"] * 100 for r in report["per_exceedance"]]
    syn_p50 = [r["syn_p50"] for r in report["per_exceedance"]]
    hist_p50 = [r["hist_p50"] for r in report["per_exceedance"]]
    syn_p10 = [r["syn_p10"] for r in report["per_exceedance"]]
    hist_p10 = [r["hist_p10"] for r in report["per_exceedance"]]

    ax.plot(xs, hist_p50, "-o", color="#1f77b4", label="Historical p50")
    ax.plot(xs, hist_p10, "--", color="#1f77b4", label="Historical p10")
    ax.plot(xs, syn_p50, "-o", color="#ff7f0e", label="Drought-subset p50")
    ax.plot(xs, syn_p10, "--", color="#ff7f0e", label="Drought-subset p10")

    # Shade fail slots
    for r in report["per_exceedance"]:
        if not r["pass"]:
            ax.axvspan(r["exceedance"] * 100 - 1, r["exceedance"] * 100 + 1,
                       alpha=0.12, color="red")

    ax.set_yscale("log")
    ax.set_xlabel("Exceedance (%)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title(f"Criterion 2 — Low-flow directionality ({report['status']})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c3_nominal(
    report: Dict,
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    status = report.get("status", "?")

    if status == "insufficient_data":
        for ax in axes:
            ax.text(0.5, 0.5, f"insufficient_data\nn={report.get('n_nominal_subset', 0)}",
                    ha="center", va="center", transform=ax.transAxes)
        fig.suptitle(f"Criterion 3 — Nominal-subset fidelity ({status})")
        fig.tight_layout()
        _save_fig(fig, output_path)
        return

    # FDC panel
    fdc = report["fdc_high_flow"]
    axes[0].bar([str(e) for e in fdc["exceedances"]],
                fdc["pct_inside_per_exceedance"], color="#2ca02c")
    axes[0].axhline(90, color="red", linestyle="--", label="90% target")
    axes[0].set_ylabel("% of nominal traces inside envelope")
    axes[0].set_title(f"High-flow FDC containment ({'pass' if fdc['pass'] else 'fail'})")
    axes[0].legend()

    # Lag-1 AC
    ac = report["lag1_ac"]
    axes[1].barh(["Synthetic", "Historical"],
                 [ac["syn_5_95"][1] - ac["syn_5_95"][0],
                  ac["hist_5_95"][1] - ac["hist_5_95"][0]],
                 left=[ac["syn_5_95"][0], ac["hist_5_95"][0]],
                 color=["#ff7f0e", "#1f77b4"], alpha=0.7)
    axes[1].set_xlabel("lag-1 AC")
    axes[1].set_title(f"Lag-1 AC 5-95 interval ({'pass' if ac['intervals_overlap'] else 'fail'})")

    # Seasonal cycle containment
    sc = report["seasonal_cycle"]
    axes[2].bar(["Monthly mean", "Monthly std"],
                [sc["min_month_pct_inside_mean"], sc["min_month_pct_inside_std"]],
                color=["#2ca02c", "#2ca02c"])
    axes[2].axhline(90, color="red", linestyle="--", label="90% target")
    axes[2].set_ylabel("Worst-month % inside envelope")
    axes[2].set_title(f"Seasonal cycle ({'pass' if sc['pass'] else 'fail'})")
    axes[2].legend()

    fig.suptitle(f"Criterion 3 — Nominal-subset fidelity ({status})")
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c4_spread(
    A_syn: np.ndarray,
    A_hist: np.ndarray,
    drought_metrics: np.ndarray,
    drought_mask: np.ndarray,
    nominal_mask: np.ndarray,
    report: Dict,
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    bins = 30

    axes[0].hist(A_hist, bins=bins, alpha=0.5, color="#1f77b4",
                 label=f"Historical blocks (std={report['std_hist']:.0f})")
    axes[0].hist(A_syn, bins=bins, alpha=0.5, color="#ff7f0e",
                 label=f"Pareto archive (std={report['std_syn']:.0f})")
    axes[0].set_xlabel("Annual mean flow (cfs-month)")
    axes[0].set_ylabel("Count")
    spread_pass_str = "pass" if report["spread_pass"] else "fail"
    axes[0].set_title(
        f"Spread ratio = {report['spread_ratio']:.2f} [0.9, 1.1] — {spread_pass_str}"
    )
    axes[0].legend()

    # Scatter annual mean vs mean_duration with corners highlighted
    axes[1].scatter(drought_metrics[:, 0], A_syn,
                    s=6, c="#888888", alpha=0.3, label="All Pareto")
    if nominal_mask.any():
        axes[1].scatter(drought_metrics[nominal_mask, 0], A_syn[nominal_mask],
                        s=12, c="#2ca02c", alpha=0.6,
                        label=f"Nominal corner (n={int(nominal_mask.sum())})")
    if drought_mask.any():
        axes[1].scatter(drought_metrics[drought_mask, 0], A_syn[drought_mask],
                        s=12, c="#d62728", alpha=0.5,
                        label=f"Drought corner (n={int(drought_mask.sum())})")
    hist_med = report["corner_bifurcation"]["hist_median_annual_mean"]
    axes[1].axhline(hist_med, color="black", linestyle="--",
                    label=f"Historical median = {hist_med:.0f}")
    axes[1].set_xlabel("mean_duration")
    axes[1].set_ylabel("Annual mean flow")
    axes[1].set_title("Corner bifurcation")
    axes[1].legend(fontsize=8)

    fig.suptitle(f"Criterion 4 — Flow-space spread ({report['status']})")
    fig.tight_layout()
    _save_fig(fig, output_path)


def plot_c5_onset(report: Dict, output_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5))
    if report.get("status") == "insufficient_data":
        ax.text(0.5, 0.5, "insufficient_data",
                ha="center", va="center", transform=ax.transAxes)
        _save_fig(fig, output_path)
        return

    months = np.arange(1, 13)
    syn = np.array(report["synthetic_histogram"], dtype=float)
    hist = np.array(report["historical_histogram"], dtype=float)
    # Normalize
    syn_n = syn / max(syn.sum(), 1.0)
    hist_n = hist / max(hist.sum(), 1.0)

    width = 0.4
    ax.bar(months - width / 2, hist_n, width, color="#1f77b4",
           label=f"Historical (n={int(hist.sum())})")
    ax.bar(months + width / 2, syn_n, width, color="#ff7f0e",
           label=f"Synthetic (n={int(syn.sum())})")

    ax.set_xticks(months)
    ax.set_xlabel("Calendar month")
    ax.set_ylabel("Fraction of drought events")
    ax.set_title(
        f"Criterion 5 — Drought onset seasonality "
        f"(χ²={report['chi2_statistic']:.1f}, p={report['p_value']:.3f} → {report['status']})"
    )
    ax.legend()
    fig.tight_layout()
    _save_fig(fig, output_path)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results_json", type=Path,
                    help="Path to results.json from a MOEA-FIND run.")
    ap.add_argument("--output-dir", type=Path, default=None,
                    help="Directory for verification_report.json and figures/. "
                         "Defaults to the results.json parent directory.")
    ap.add_argument("--min-subset-size", type=int, default=100,
                    help="Minimum traces per subset for Criteria 2, 3, 5 "
                         "(default: 100).")
    args = ap.parse_args()

    if not args.results_json.exists():
        raise FileNotFoundError(args.results_json)
    variant_dir = args.results_json.parent
    output_dir = args.output_dir or variant_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"

    print(f"[verify] loading {args.results_json}")
    results = json.loads(args.results_json.read_text())
    pareto_traces_1d = [np.asarray(t, dtype=float) for t in results["pareto_traces_1d"]]
    pareto_traces_2d = [np.asarray(t, dtype=float) for t in results["pareto_traces_2d"]]
    drought_metrics = np.asarray(results["drought_metrics"], dtype=float)
    objective_keys = tuple(results.get(
        "objective_keys", ["mean_duration", "mean_avg_severity"]))
    anti_ideal = np.asarray(results.get("anti_ideal", []), dtype=float)
    n_years = int(results["n_years_out"])
    ssi_timescale = int(results.get("ssi_timescale", 3))

    print(f"[verify] {len(pareto_traces_1d)} Pareto traces, "
          f"K={drought_metrics.shape[1]}, T={n_years} years, SSI-{ssi_timescale}")

    # Historical data + blocks (match script 04 conventions) ------------------
    print("[verify] loading historical flows + building T-year blocks")
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir=cache_dir)
    hist_blocks_1d = resample_historical_blocks(monthly_1d, T_years=n_years, stride=1)
    hist_blocks_2d = resample_historical_blocks_2d(monthly_2d, T_years=n_years, stride=1)
    print(f"[verify] {len(hist_blocks_1d)} historical blocks at T={n_years}")

    # Historical block drought chars — load cache if present, else recompute --
    hist_char_npz = variant_dir / "historical_block_chars.npz"
    ssi_calc = None
    if hist_char_npz.exists():
        print(f"[verify] loading cached {hist_char_npz}")
        npz = np.load(hist_char_npz, allow_pickle=True)
        hist_block_chars = np.asarray(npz["chars"], dtype=float)
        cached_keys = [str(k) for k in npz["objective_keys"]]
        if tuple(cached_keys) != objective_keys:
            warnings.warn(
                f"cached objective_keys {cached_keys} != results "
                f"{objective_keys}; recomputing"
            )
            hist_block_chars = None

    if ssi_calc is None:
        # Fit the same SSI calculator the MOEA used — prefit on full history.
        _ssi_series, ssi_calc = compute_ssi(monthly_1d, timescale=ssi_timescale)

    if not hist_char_npz.exists() or hist_block_chars is None:
        print("[verify] recomputing historical block drought chars")
        hist_block_chars = compute_historical_block_chars(
            monthly_1d, T_years=n_years, ssi_calc=ssi_calc,
            objective_keys=objective_keys, stride=1,
        )

    # Subsets ----------------------------------------------------------------
    hist_medians = np.median(hist_block_chars, axis=0)
    drought_mask = (drought_metrics >= hist_medians).any(axis=1)
    nominal_mask = (drought_metrics < hist_medians).all(axis=1)
    n_drought = int(drought_mask.sum())
    n_nominal = int(nominal_mask.sum())
    print(f"[verify] subsets: drought={n_drought}, nominal={n_nominal} "
          f"(union vs intersection complement on hist medians={hist_medians})")

    drought_traces_1d = [pareto_traces_1d[i]
                         for i in np.where(drought_mask)[0]]
    nominal_traces_1d = [pareto_traces_1d[i]
                         for i in np.where(nominal_mask)[0]]
    nominal_traces_2d = [pareto_traces_2d[i]
                         for i in np.where(nominal_mask)[0]]

    # Criteria ---------------------------------------------------------------
    print("[verify] criterion 1 — drought-space coverage")
    c1 = criterion_1_coverage(drought_metrics, hist_block_chars, objective_keys)

    print("[verify] criterion 2 — low-flow directionality")
    c2 = criterion_2_low_flow_directionality(drought_traces_1d, hist_blocks_1d)

    print("[verify] criterion 3 — non-drought fidelity")
    c3 = criterion_3_nominal_fidelity(
        nominal_traces_1d, nominal_traces_2d, hist_blocks_1d, hist_blocks_2d,
    )

    print("[verify] criterion 4 — flow-space spread")
    c4 = criterion_4_flow_space_spread(
        pareto_traces_1d, hist_blocks_1d, drought_metrics,
        hist_block_chars, drought_mask, nominal_mask,
    )

    print("[verify] criterion 5 — drought onset seasonality (SSI re-transform)")
    c5 = criterion_5_drought_onset(
        drought_traces_1d, hist_blocks_1d, ssi_calc,
        start_date="2100-01-01",
    )

    # Assemble report --------------------------------------------------------
    report = {
        "archive": str(args.results_json),
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

    report_path = output_dir / "verification_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"[verify] wrote {report_path} — overall={report['overall_status']}")

    # Figures ----------------------------------------------------------------
    plot_c1_coverage(drought_metrics, hist_block_chars,
                     objective_keys, anti_ideal,
                     fig_dir / "fig_verify_c1_drought_coverage.pdf")
    plot_c2_low_flow(c2, fig_dir / "fig_verify_c2_drought_subset_fdc.pdf")
    plot_c3_nominal(c3, fig_dir / "fig_verify_c3_nominal_subset.pdf")
    # Prepare annual-mean arrays for c4 figure
    A_syn = np.array([_annual_totals(t).mean() for t in pareto_traces_1d])
    A_hist = np.array([_annual_totals(b).mean() for b in hist_blocks_1d])
    plot_c4_spread(A_syn, A_hist, drought_metrics,
                   drought_mask, nominal_mask, c4,
                   fig_dir / "fig_verify_c4_annual_spread.pdf")
    plot_c5_onset(c5, fig_dir / "fig_verify_c5_drought_onset.pdf")

    # Compact terminal summary -----------------------------------------------
    print("\n[verify] summary")
    for key, crit in report["criteria"].items():
        print(f"  {key:25s} {crit['status']}")
    print(f"  {'overall':25s} {report['overall_status']}")


if __name__ == "__main__":
    main()
