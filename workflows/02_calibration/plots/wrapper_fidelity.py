"""Render Kirsch-wrapper fidelity figures (SI C).

Reads ``outputs/02_calibration/wrapper_fidelity/{ensembles_2d.npz,
drought_clouds.npz}`` and writes the Phase A (trace-level) and Phase B
(drought-space) PDFs into
``figures/02_calibration/wrapper_fidelity/``.

Plotting-only -- never re-runs trace generation.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.style import apply_style, WATER_YEAR_MONTHS  # noqa: E402

STAGE = "02_calibration"
DRIVER = "wrapper_fidelity"

ENS_COLORS: Dict[str, str] = {
    "historical": "#000000",
    "kirsch_baseline": "#7f7f7f",
    "index": "#2b6cb0",
    "residual": "#c05621",
}
ENS_LABELS: Dict[str, str] = {
    "historical": "historical T-blocks",
    "kirsch_baseline": "Kirsch baseline (no DV)",
    "index": "index wrapper (U[0,1] DVs)",
    "residual": "residual wrapper (U[0,1] DVs)",
}
ENS_ORDER: Tuple[str, ...] = ("historical", "kirsch_baseline", "index", "residual")


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def _acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    n = len(x)
    mean = np.mean(x)
    var = np.var(x)
    if var < 1e-12:
        return np.zeros(max_lag + 1)
    out = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        out[lag] = np.mean((x[: n - lag] - mean) * (x[lag:] - mean)) / var
    return out


def _envelope(stack: np.ndarray, pcts=(10, 50, 90)):
    return tuple(np.percentile(stack, p, axis=0) for p in pcts)


def _stack_fdc(traces_1d: List[np.ndarray], grid: np.ndarray) -> np.ndarray:
    stacked = np.empty((len(traces_1d), len(grid)))
    for i, tr in enumerate(traces_1d):
        srt = np.sort(tr)[::-1]
        ex = np.arange(1, len(srt) + 1) / (len(srt) + 1)
        stacked[i] = np.interp(grid, ex, srt)
    return stacked


def _cross_month_corr(traces_2d: np.ndarray) -> np.ndarray:
    mats = [np.corrcoef(t, rowvar=False) for t in traces_2d if t.shape[0] >= 2]
    return np.mean(np.array(mats), axis=0) if mats else np.zeros((12, 12))


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------

def _fig_baseline_check(ens_1d, ens_2d, fig_path):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    ax_fdc, ax_lag1, ax_ann = axes

    any_trace = next(iter(ens_1d.values()))[0]
    grid = np.arange(1, len(any_trace) + 1) / (len(any_trace) + 1)
    grid_pct = grid * 100
    for ens in ENS_ORDER:
        traces = ens_1d[ens]
        if not traces:
            continue
        stack = _stack_fdc(traces, grid)
        lo, mid, hi = _envelope(stack)
        ax_fdc.fill_between(grid_pct, lo, hi, alpha=0.18, color=ENS_COLORS[ens],
                            label=f"{ENS_LABELS[ens]} (n={len(traces)})")
        ax_fdc.semilogy(grid_pct, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax_fdc.set_xlabel("Exceedance probability (%)")
    ax_fdc.set_ylabel("Flow (cfs, log)")
    ax_fdc.set_title("(a) FDC 10-90% envelope")
    ax_fdc.legend(fontsize=7, loc="upper right", framealpha=0.85)

    for ens in ENS_ORDER:
        traces = ens_1d[ens]
        if not traces:
            continue
        lag1 = np.array([_acf(t, 1)[1] for t in traces])
        ax_lag1.hist(lag1, bins=40, density=True, alpha=0.45,
                     color=ENS_COLORS[ens],
                     label=f"{ENS_LABELS[ens]} (med={np.median(lag1):.2f})")
    ax_lag1.set_xlabel("Lag-1 autocorrelation")
    ax_lag1.set_ylabel("Density")
    ax_lag1.set_title("(b) Lag-1 ACF distribution")
    ax_lag1.legend(fontsize=7, framealpha=0.85)

    for ens in ENS_ORDER:
        traces = ens_2d[ens]
        if traces.shape[0] == 0:
            continue
        totals = traces.sum(axis=2).flatten()
        ax_ann.hist(totals, bins=60, density=True, alpha=0.4,
                    color=ENS_COLORS[ens], label=ENS_LABELS[ens])
    ax_ann.set_xlabel("Annual total flow (cfs)")
    ax_ann.set_ylabel("Density")
    ax_ann.set_title("(c) Annual total distribution")
    ax_ann.legend(fontsize=7, framealpha=0.85)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_monthly_marginal_boxplots(ens_2d, fig_path):
    fig, ax = plt.subplots(figsize=(12.0, 4.5))
    months = np.arange(12)
    width = 0.2
    active = [e for e in ENS_ORDER if ens_2d[e].shape[0] > 0]
    n_active = len(active)
    start = -(n_active - 1) / 2.0

    for i, ens in enumerate(active):
        data = ens_2d[ens].reshape(-1, 12)  # (n_traces*n_years, 12)
        positions = months + (start + i) * width
        box_data = [data[:, m] for m in range(12)]
        bp = ax.boxplot(box_data, positions=positions, widths=width * 0.9,
                        patch_artist=True, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(ENS_COLORS[ens])
            patch.set_alpha(0.45)
            patch.set_edgecolor(ENS_COLORS[ens])
        for w in bp["whiskers"]:
            w.set_color(ENS_COLORS[ens])
        for c in bp["caps"]:
            c.set_color(ENS_COLORS[ens])
        for m in bp["medians"]:
            m.set_color("black")
            m.set_linewidth(1.2)

    ax.set_xticks(months)
    ax.set_xticklabels(WATER_YEAR_MONTHS)
    ax.set_ylabel("Monthly flow (cfs)")
    ax.set_yscale("log")
    ax.set_title("Per-month marginals by ensemble (log scale)")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ENS_COLORS[e], alpha=0.45,
                              edgecolor=ENS_COLORS[e], label=ENS_LABELS[e])
               for e in active]
    ax.legend(handles=handles, fontsize=8, loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_fdc_envelopes(ens_1d, fig_path):
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    any_trace = next(iter(ens_1d.values()))[0]
    grid = np.arange(1, len(any_trace) + 1) / (len(any_trace) + 1)
    grid_pct = grid * 100
    for ens in ENS_ORDER:
        traces = ens_1d[ens]
        if not traces:
            continue
        stack = _stack_fdc(traces, grid)
        lo, mid, hi = _envelope(stack)
        ax.fill_between(grid_pct, lo, hi, alpha=0.2, color=ENS_COLORS[ens],
                        label=f"{ENS_LABELS[ens]} 10-90% (n={len(traces)})")
        ax.semilogy(grid_pct, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax.set_xlabel("Exceedance probability (%)")
    ax.set_ylabel("Flow (cfs, log)")
    ax.set_title("Flow duration curve envelopes")
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_acf_envelopes(ens_1d, fig_path, max_lag=24):
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    lags = np.arange(max_lag + 1)
    per_ens: Dict[str, np.ndarray] = {}
    for ens in ENS_ORDER:
        traces = ens_1d[ens]
        if not traces:
            continue
        acfs = np.array([_acf(t, max_lag) for t in traces])
        per_ens[ens] = acfs
        lo, mid, hi = _envelope(acfs)
        ax.fill_between(lags, lo, hi, alpha=0.18, color=ENS_COLORS[ens],
                        label=f"{ENS_LABELS[ens]} 10-90%")
        ax.plot(lags, mid, color=ENS_COLORS[ens], linewidth=1.6)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Lag (months)")
    ax.set_ylabel("Autocorrelation")
    ax.set_title(f"ACF envelopes (lags 0-{max_lag})")
    ax.legend(fontsize=8, framealpha=0.9)

    ins = ax.inset_axes([0.55, 0.50, 0.42, 0.42])
    positions = np.arange(2)
    width = 0.18
    active = [e for e in ENS_ORDER if e in per_ens]
    for i, ens in enumerate(active):
        arr = per_ens[ens]
        medians = [np.median(arr[:, 1]), np.median(arr[:, 12])]
        lows = [np.percentile(arr[:, 1], 10), np.percentile(arr[:, 12], 10)]
        his = [np.percentile(arr[:, 1], 90), np.percentile(arr[:, 12], 90)]
        offsets = (i - (len(active) - 1) / 2.0) * width
        ins.bar(positions + offsets, medians, width=width * 0.95,
                color=ENS_COLORS[ens], alpha=0.8)
        ins.errorbar(
            positions + offsets, medians,
            yerr=[np.array(medians) - np.array(lows),
                  np.array(his) - np.array(medians)],
            fmt="none", ecolor="black", linewidth=0.8, capsize=2,
        )
    ins.set_xticks(positions)
    ins.set_xticklabels(["lag 1", "lag 12"], fontsize=7)
    ins.tick_params(axis="y", labelsize=7)
    ins.set_title("median +/- 10/90%", fontsize=7)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_cross_month_corr_heatmap(ens_2d, fig_path):
    active = [e for e in ENS_ORDER if ens_2d[e].shape[0] > 0]
    mats = {e: _cross_month_corr(ens_2d[e]) for e in active}

    n = len(active)
    fig, axes = plt.subplots(2, n, figsize=(3.0 * n + 1.0, 6.3))
    if n == 1:
        axes = axes.reshape(2, 1)

    for j, ens in enumerate(active):
        ax = axes[0, j]
        im = ax.imshow(mats[ens], vmin=-0.5, vmax=1.0, cmap="RdBu_r", aspect="equal")
        ax.set_title(ENS_LABELS[ens], fontsize=9)
        ax.set_xticks(range(12)); ax.set_xticklabels(WATER_YEAR_MONTHS, fontsize=6, rotation=90)
        ax.set_yticks(range(12)); ax.set_yticklabels(WATER_YEAR_MONTHS, fontsize=6)
    cbar = fig.colorbar(im, ax=axes[0, :].tolist(), shrink=0.75)
    cbar.set_label("cross-month Pearson r", fontsize=8)

    ref_key = "historical" if "historical" in mats else active[0]
    ref = mats[ref_key]
    for j, ens in enumerate(active):
        ax = axes[1, j]
        diff = mats[ens] - ref
        im2 = ax.imshow(diff, vmin=-0.4, vmax=0.4, cmap="RdBu_r", aspect="equal")
        frob = np.linalg.norm(diff, ord="fro")
        ax.set_title(f"{ens} - {ref_key}  (||.||_F = {frob:.2f})", fontsize=8)
        ax.set_xticks(range(12)); ax.set_xticklabels(WATER_YEAR_MONTHS, fontsize=6, rotation=90)
        ax.set_yticks(range(12)); ax.set_yticklabels(WATER_YEAR_MONTHS, fontsize=6)
    cbar2 = fig.colorbar(im2, ax=axes[1, :].tolist(), shrink=0.75)
    cbar2.set_label("delta cross-month r", fontsize=8)

    fig.suptitle("Cross-month correlation (ensemble-averaged)", fontsize=11)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_seasonal_cycle(ens_2d, fig_path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    months = np.arange(12)
    for ens in ENS_ORDER:
        traces = ens_2d[ens]
        if traces.shape[0] == 0:
            continue
        means = traces.mean(axis=1)
        stds = traces.std(axis=1, ddof=1)
        for ax, stat in [(ax1, means), (ax2, stds)]:
            lo, mid, hi = _envelope(stat)
            ax.fill_between(months, lo, hi, alpha=0.18, color=ENS_COLORS[ens])
            ax.plot(months, mid, color=ENS_COLORS[ens], linewidth=1.6,
                    label=ENS_LABELS[ens])
    for ax, ttl, ylbl in [(ax1, "Seasonal mean", "Mean monthly flow (cfs)"),
                          (ax2, "Seasonal std", "Std monthly flow (cfs)")]:
        ax.set_xticks(months); ax.set_xticklabels(WATER_YEAR_MONTHS, rotation=45)
        ax.set_ylabel(ylbl); ax.set_title(ttl)
    ax1.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_annual_tails(ens_2d, fig_path):
    fig, (ax_tot, ax_min) = plt.subplots(1, 2, figsize=(11.0, 4.5))
    for ens in ENS_ORDER:
        traces = ens_2d[ens]
        if traces.shape[0] == 0:
            continue
        totals = traces.sum(axis=2).flatten()
        mins = traces.min(axis=2).flatten()
        for ax, data in [(ax_tot, totals), (ax_min, mins)]:
            srt = np.sort(data)
            cdf = np.arange(1, len(srt) + 1) / (len(srt) + 1)
            ax.plot(srt, cdf, color=ENS_COLORS[ens], linewidth=1.6,
                    label=ENS_LABELS[ens])
    ax_tot.set_xlabel("Annual total flow (cfs)")
    ax_tot.set_title("Annual total ECDF")
    ax_tot.set_ylabel("Cumulative probability")
    ax_min.set_xlabel("Annual minimum monthly flow (cfs)")
    ax_min.set_title("Annual min-month ECDF (drought tail)")
    ax_min.set_xscale("log")
    for ax in (ax_tot, ax_min):
        ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_drought_space_2d(mode_clouds, hist_cloud, baseline_cloud,
                          exp05_cloud, objective_labels, fig_path):
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    if exp05_cloud is not None and exp05_cloud.shape[0] > 0:
        ax.scatter(exp05_cloud[:, 0], exp05_cloud[:, 1],
                   color="#d9d9d9", s=5, alpha=0.35,
                   label=f"Kirsch library (n={len(exp05_cloud)})")
    if hist_cloud is not None and hist_cloud.shape[0] > 0:
        ax.scatter(hist_cloud[:, 0], hist_cloud[:, 1],
                   color=ENS_COLORS["historical"], s=30, marker="x",
                   label=f"historical T-blocks (n={len(hist_cloud)})")
    if baseline_cloud is not None and baseline_cloud.shape[0] > 0:
        ax.scatter(baseline_cloud[:, 0], baseline_cloud[:, 1],
                   color=ENS_COLORS["kirsch_baseline"], s=8, alpha=0.45,
                   label=f"Kirsch baseline (n={len(baseline_cloud)})")
    for mode in ("index", "residual"):
        cloud = mode_clouds.get(mode)
        if cloud is None or cloud.shape[0] == 0:
            continue
        ax.scatter(cloud[:, 0], cloud[:, 1],
                   color=ENS_COLORS[mode], s=10, alpha=0.55, marker="o",
                   edgecolor="none",
                   label=f"{ENS_LABELS[mode]} (n={len(cloud)})")
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_title("Drought-space coverage under random DVs")
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_drought_space_3d(mode_clouds, hist_cloud, baseline_cloud,
                          exp05_cloud, objective_labels, fig_path):
    fig = plt.figure(figsize=(8.0, 6.5))
    ax = fig.add_subplot(111, projection="3d")
    if exp05_cloud is not None and exp05_cloud.shape[0] > 0 and exp05_cloud.shape[1] >= 3:
        ax.scatter(exp05_cloud[:, 0], exp05_cloud[:, 1], exp05_cloud[:, 2],
                   color="#d9d9d9", s=3, alpha=0.3,
                   label=f"Kirsch library (n={len(exp05_cloud)})")
    if hist_cloud is not None and hist_cloud.shape[0] > 0:
        ax.scatter(hist_cloud[:, 0], hist_cloud[:, 1], hist_cloud[:, 2],
                   color=ENS_COLORS["historical"], s=20, marker="x",
                   label=f"historical T-blocks (n={len(hist_cloud)})")
    if baseline_cloud is not None and baseline_cloud.shape[0] > 0:
        ax.scatter(baseline_cloud[:, 0], baseline_cloud[:, 1], baseline_cloud[:, 2],
                   color=ENS_COLORS["kirsch_baseline"], s=6, alpha=0.4,
                   label=f"Kirsch baseline (n={len(baseline_cloud)})")
    for mode in ("index", "residual"):
        cloud = mode_clouds.get(mode)
        if cloud is None or cloud.shape[0] == 0 or cloud.shape[1] < 3:
            continue
        ax.scatter(cloud[:, 0], cloud[:, 1], cloud[:, 2],
                   color=ENS_COLORS[mode], s=8, alpha=0.5,
                   label=f"{ENS_LABELS[mode]} (n={len(cloud)})")
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_zlabel(objective_labels[2])
    ax.set_title("Drought-space 3D coverage (random DVs)")
    ax.legend(fontsize=7, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig_coverage_vs_n(mode_clouds, exp05_cloud, objective_labels, fig_path,
                       ns=(50, 100, 200, 500, 1000)):
    n_obj = len(objective_labels)
    fig, axes = plt.subplots(1, n_obj, figsize=(4.0 * n_obj, 4.2), sharex=True)
    if n_obj == 1:
        axes = [axes]
    for j, lbl in enumerate(objective_labels):
        ax = axes[j]
        for mode in ("index", "residual"):
            cloud = mode_clouds.get(mode)
            if cloud is None or cloud.shape[0] == 0 or cloud.shape[1] <= j:
                continue
            max_curve = []
            for n in ns:
                n_use = min(n, cloud.shape[0])
                max_curve.append(np.max(cloud[:n_use, j]))
            ax.plot(ns, max_curve, "o-", color=ENS_COLORS[mode],
                    label=ENS_LABELS[mode], linewidth=1.6, markersize=5)
        if exp05_cloud is not None and exp05_cloud.shape[0] > 0 and exp05_cloud.shape[1] > j:
            ax.axhline(np.max(exp05_cloud[:, j]), color="#7f7f7f",
                       linestyle="--", linewidth=1.2,
                       label=f"library max (n={len(exp05_cloud)})")
        ax.set_xlabel("Ensemble size N")
        ax.set_ylabel(f"max {lbl}")
        ax.set_xscale("log")
        ax.set_title(lbl)
        ax.legend(fontsize=7, framealpha=0.9)
    fig.suptitle("Per-axis range coverage vs random-DV ensemble size", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    apply_style()

    in_dir = stage_output_dir(STAGE, DRIVER, create=False)
    ens_path = in_dir / "ensembles_2d.npz"
    if not ens_path.exists():
        sys.exit(f"missing {ens_path} -- run the compute driver first")
    npz = np.load(ens_path)
    ens_2d = {k: npz[k] for k in ENS_ORDER}
    ens_1d = {k: [arr.flatten() for arr in ens_2d[k]] for k in ENS_ORDER}

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    print(f"[plots/02/wrapper_fidelity] fig_dir={fig_dir}")
    _fig_baseline_check(ens_1d, ens_2d, fig_dir / "fig_wrapper_baseline_check.pdf")
    _fig_monthly_marginal_boxplots(ens_2d, fig_dir / "fig_monthly_marginal_boxplots.pdf")
    _fig_fdc_envelopes(ens_1d, fig_dir / "fig_fdc_envelopes.pdf")
    _fig_acf_envelopes(ens_1d, fig_dir / "fig_acf_envelopes.pdf")
    _fig_cross_month_corr_heatmap(ens_2d, fig_dir / "fig_cross_month_corr_heatmap.pdf")
    _fig_seasonal_cycle(ens_2d, fig_dir / "fig_seasonal_cycle.pdf")
    _fig_annual_tails(ens_2d, fig_dir / "fig_annual_tails.pdf")

    clouds_path = in_dir / "drought_clouds.npz"
    if not clouds_path.exists():
        print(f"  drought_clouds.npz absent -- skipping Phase B figures.")
        return
    cnpz = np.load(clouds_path, allow_pickle=True)
    mode_clouds = {"index": cnpz["index"], "residual": cnpz["residual"]}
    hist_cloud = cnpz["historical"]
    baseline_cloud = cnpz["kirsch_baseline"]
    exp05_cloud = cnpz["exp05"]
    objective_labels = tuple(str(s) for s in cnpz["objective_labels"])

    _fig_drought_space_2d(
        {m: c[:, :2] for m, c in mode_clouds.items()},
        hist_cloud[:, :2] if hist_cloud.shape[0] > 0 else hist_cloud,
        baseline_cloud[:, :2],
        exp05_cloud[:, :2] if exp05_cloud.shape[0] > 0 else exp05_cloud,
        objective_labels[:2],
        fig_dir / "fig_drought_space_random_2d.pdf",
    )
    if len(objective_labels) >= 3:
        _fig_drought_space_3d(
            mode_clouds, hist_cloud, baseline_cloud, exp05_cloud,
            objective_labels, fig_dir / "fig_drought_space_random_3d.pdf",
        )
    _fig_coverage_vs_n(
        mode_clouds, exp05_cloud, objective_labels,
        fig_dir / "fig_coverage_vs_ensemble_size.pdf",
    )
    print(f"  done.")


if __name__ == "__main__":
    main()
