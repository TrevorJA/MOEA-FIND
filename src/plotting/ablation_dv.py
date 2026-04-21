"""Plotting composition helpers for the DV-uniformity ablation SI figures.

Wraps the existing primitives in :mod:`src.plotting.trace_diagnostics` so the
compare script (``workflows/experiments/14_dv_uniformity_compare.py``) can
render side-by-side or overlaid comparisons of the hydrologic arm and the
DV-uniform arm against the same historical block envelope, without
duplicating the envelope-percentile machinery.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from src.plotting.style import COLORS, WATER_YEAR_MONTHS, apply_style
from src.plotting.trace_diagnostics import _acf, _envelope, _fdc


ARM_COLORS: Dict[str, str] = {
    "hydrologic": "#2b6cb0",  # steel blue
    "dv_uniform": "#c05621",  # terracotta (legacy aggregate)
    "dv_l2_star": "#c05621",  # terracotta
    "dv_ks": "#6b46c1",       # purple
    "dv_ad": "#2f855a",       # forest green
    # wrapper-mode ablation (exp15/16)
    "index":    "#2b6cb0",    # steel blue
    "residual": "#c05621",    # terracotta
}

ARM_LABELS: Dict[str, str] = {
    "hydrologic": "hydrologic (5 constraints)",
    "dv_uniform": "DV uniform",
    "dv_l2_star": "DV L2-star",
    "dv_ks": "DV KS",
    "dv_ad": "DV Anderson-Darling",
    # wrapper-mode ablation (exp15/16)
    "index":    "Index wrapper",
    "residual": "Residual wrapper",
}


def _stack_fdc(traces: Sequence[np.ndarray], grid: np.ndarray) -> np.ndarray:
    stacked = np.empty((len(traces), len(grid)))
    for i, tr in enumerate(traces):
        e, f = _fdc(tr)
        stacked[i] = np.interp(grid, e, f)
    return stacked


def plot_hydrology_panels_two_arms(
    traces_1d_by_arm: Dict[str, List[np.ndarray]],
    traces_2d_by_arm: Dict[str, List[np.ndarray]],
    historical_blocks_1d: List[np.ndarray],
    historical_blocks_2d: List[np.ndarray],
    max_lag: int = 24,
    figsize: Tuple[float, float] = (11.0, 7.0),
) -> Tuple[plt.Figure, np.ndarray]:
    """Four-panel ACF / FDC / seasonal-mean / seasonal-std comparison.

    Draws the historical block envelope once, then overlays one envelope
    per arm keyed by ``ARM_COLORS``. Mirrors
    :func:`src.plotting.trace_diagnostics.plot_hydrology_panels` but for
    two synthetic ensembles on the same axes.
    """
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    ax_acf, ax_fdc = axes[0]
    ax_mean, ax_std = axes[1]

    # -- (a) ACF envelope --
    lags = np.arange(max_lag + 1)
    hist_acfs = np.array([_acf(b, max_lag) for b in historical_blocks_1d])
    h_lo, h_mid, h_hi = _envelope(hist_acfs)
    ax_acf.fill_between(lags, h_lo, h_hi, alpha=0.25,
                        color=COLORS["historical"],
                        label=f"historical blocks (n={len(historical_blocks_1d)}, 10-90%)")
    ax_acf.plot(lags, h_mid, color=COLORS["historical"],
                linewidth=2, label="historical median", zorder=10)
    for arm, traces in traces_1d_by_arm.items():
        if not traces:
            continue
        acfs = np.array([_acf(t, max_lag) for t in traces])
        s_lo, s_mid, s_hi = _envelope(acfs)
        ax_acf.fill_between(lags, s_lo, s_hi, alpha=0.18,
                            color=ARM_COLORS.get(arm, "gray"),
                            label=f"{ARM_LABELS.get(arm, arm)} (n={len(traces)}, 10-90%)")
        ax_acf.plot(lags, s_mid, "--", color=ARM_COLORS.get(arm, "gray"),
                    linewidth=1.5, label=f"{ARM_LABELS.get(arm, arm)} median")
    ax_acf.axhline(0, color="gray", linewidth=0.5, linestyle=":")
    ax_acf.set_xlabel("lag (months)")
    ax_acf.set_ylabel("autocorrelation")
    ax_acf.set_title("(a) autocorrelation")
    ax_acf.legend(fontsize=7, loc="best", framealpha=0.9)

    # -- (b) FDC envelope (log-y) --
    # Grid sized to the shortest trace (defensive).
    ref_traces = next(iter(traces_1d_by_arm.values()), historical_blocks_1d)
    n_months = len(ref_traces[0]) if ref_traces else len(historical_blocks_1d[0])
    grid = np.arange(1, n_months + 1) / (n_months + 1)
    grid_pct = grid * 100
    hist_stack = _stack_fdc(historical_blocks_1d, grid)
    h_lo, h_mid, h_hi = _envelope(hist_stack)
    ax_fdc.fill_between(grid_pct, h_lo, h_hi, alpha=0.25,
                        color=COLORS["historical"],
                        label=f"historical (n={len(historical_blocks_1d)}, 10-90%)")
    ax_fdc.semilogy(grid_pct, h_mid, color=COLORS["historical"],
                    linewidth=2, label="historical median", zorder=10)
    for arm, traces in traces_1d_by_arm.items():
        if not traces:
            continue
        stack = _stack_fdc(traces, grid)
        s_lo, s_mid, s_hi = _envelope(stack)
        ax_fdc.fill_between(grid_pct, s_lo, s_hi, alpha=0.18,
                            color=ARM_COLORS.get(arm, "gray"),
                            label=f"{ARM_LABELS.get(arm, arm)} (n={len(traces)}, 10-90%)")
        ax_fdc.semilogy(grid_pct, s_mid, "--",
                        color=ARM_COLORS.get(arm, "gray"),
                        linewidth=1.5, label=f"{ARM_LABELS.get(arm, arm)} median")
    ax_fdc.set_xlabel("exceedance probability (%)")
    ax_fdc.set_ylabel("flow (cfs)")
    ax_fdc.set_title("(b) flow duration curve")
    ax_fdc.legend(fontsize=7, loc="best", framealpha=0.9)

    # -- (c) Seasonal mean; (d) seasonal std --
    months = np.arange(12)
    h_means = np.array([b.mean(axis=0) for b in historical_blocks_2d])
    h_stds = np.array([b.std(axis=0, ddof=1) for b in historical_blocks_2d])
    for ax, hist_stat, ylabel, title, key in [
        (ax_mean, h_means, "mean monthly flow (cfs)",
         "(c) seasonal mean", "mean"),
        (ax_std, h_stds, "std monthly flow (cfs)",
         "(d) seasonal variability", "std"),
    ]:
        p_lo, p_mid, p_hi = _envelope(hist_stat)
        ax.fill_between(months, p_lo, p_hi, alpha=0.25,
                        color=COLORS["historical"],
                        label=f"historical (n={len(hist_stat)}, 10-90%)")
        ax.plot(months, p_mid, "o-", color=COLORS["historical"],
                linewidth=2, markersize=5,
                label="historical median", zorder=10)
        for arm, traces in traces_2d_by_arm.items():
            if not traces:
                continue
            if key == "mean":
                syn_stat = np.array([t.mean(axis=0) for t in traces])
            else:
                syn_stat = np.array([t.std(axis=0, ddof=1) for t in traces])
            s_lo, s_mid, s_hi = _envelope(syn_stat)
            ax.fill_between(months, s_lo, s_hi, alpha=0.18,
                            color=ARM_COLORS.get(arm, "gray"),
                            label=f"{arm} (n={len(syn_stat)}, 10-90%)")
            ax.plot(months, s_mid, "--", color=ARM_COLORS.get(arm, "gray"),
                    linewidth=1.5, label=f"{ARM_LABELS.get(arm, arm)} median")
        ax.set_xticks(months)
        ax.set_xticklabels(WATER_YEAR_MONTHS, rotation=45, fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best", framealpha=0.9)

    fig.tight_layout()
    return fig, axes


def plot_pareto_2d_overlay(
    pareto_by_arm: Dict[str, np.ndarray],
    historical_cloud: np.ndarray,
    historical_point: Tuple[float, float],
    objective_labels: Tuple[str, str] = (
        "mean duration (months)", "mean avg severity"),
    figsize: Tuple[float, float] = (7.0, 6.0),
) -> Tuple[plt.Figure, plt.Axes]:
    """Single-panel 2D Pareto overlay, one colour per arm (legacy)."""
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)
    if historical_cloud is not None and len(historical_cloud) > 0:
        ax.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                   s=15, alpha=0.25, color="gray",
                   label=f"historical blocks (n={len(historical_cloud)})",
                   zorder=1)
    for arm, dm in pareto_by_arm.items():
        if dm is None or len(dm) == 0:
            continue
        ax.scatter(dm[:, 0], dm[:, 1],
                   s=18, alpha=0.45,
                   color=ARM_COLORS.get(arm, "gray"),
                   label=f"{ARM_LABELS.get(arm, arm)} (n={len(dm)})",
                   zorder=3)
    ax.scatter(*historical_point, marker="*", s=260, color="black",
               zorder=6, label="historical point")
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_title("Pareto drought-space overlay across constraint formulations")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    return fig, ax


def plot_pareto_2d_subpanels(
    pareto_by_arm: Dict[str, np.ndarray],
    historical_cloud: np.ndarray,
    historical_median: Tuple[float, float],
    kirsch_cloud: Optional[np.ndarray] = None,
    objective_labels: Tuple[str, str] = (
        "mean duration (months)", "mean avg severity"),
    arm_order: Optional[List[str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
) -> Tuple[plt.Figure, np.ndarray]:
    """Multi-panel 2D Pareto, one panel per arm plus a reference panel.

    All panels share the same x/y axis limits and gridlines for easy
    cross-arm comparison. The first panel (reference) shows the Kirsch
    library cloud and historical T-blocks. Each arm panel shows the arm's
    Pareto front with the historical T-block cloud as a faint backdrop.

    ``historical_median`` should be the element-wise median of
    ``historical_cloud`` (a T-block scale point, not the full-record
    characteristics which span a different number of drought events).
    """
    apply_style()
    arms = arm_order if arm_order is not None else list(pareto_by_arm.keys())
    n_panels = len(arms) + 1  # +1 for reference

    ncols = min(3, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    if figsize is None:
        figsize = (4.5 * ncols, 4.2 * nrows)

    fig, axes_2d = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = axes_2d.flatten()

    # Global axis limits across all data sources
    all_x, all_y = [], []
    for src in [historical_cloud, kirsch_cloud]:
        if src is not None and len(src) > 0:
            all_x.extend(src[:, 0])
            all_y.extend(src[:, 1])
    for dm in pareto_by_arm.values():
        if dm is not None and len(dm) > 0:
            all_x.extend(dm[:, 0])
            all_y.extend(dm[:, 1])
    if not all_x:
        return fig, axes_flat
    xpad = max(0.05 * (max(all_x) - min(all_x)), 0.1)
    ypad = max(0.05 * (max(all_y) - min(all_y)), 0.05)
    xlim = (min(all_x) - xpad, max(all_x) + xpad)
    ylim = (min(all_y) - ypad, max(all_y) + ypad)

    panel_letters = "abcdefghijklmnop"
    n_hist = len(historical_cloud) if historical_cloud is not None else 0
    n_kirsch = len(kirsch_cloud) if kirsch_cloud is not None else 0

    def _setup_ax(ax, title):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_xlabel(objective_labels[0], fontsize=9)
        ax.set_ylabel(objective_labels[1], fontsize=9)
        ax.set_title(title, fontsize=10)

    # --- Panel 0: Reference ---
    ax0 = axes_flat[0]
    if kirsch_cloud is not None and n_kirsch > 0:
        kc = kirsch_cloud
        if n_kirsch > 2000:
            rng = np.random.default_rng(0)
            kc = kc[rng.choice(n_kirsch, size=2000, replace=False)]
        ax0.scatter(kc[:, 0], kc[:, 1], s=6, alpha=0.12, color="#a0c4e8", zorder=1)
    if historical_cloud is not None and n_hist > 0:
        ax0.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                    s=8, color="black", zorder=10)
    ax0.scatter(*historical_median, marker="*", s=160, color="black", zorder=11)
    _setup_ax(ax0, f"({panel_letters[0]}) Kirsch reference + historical T-blocks")

    # --- Arm panels ---
    for i, arm in enumerate(arms):
        ax = axes_flat[i + 1]
        dm = pareto_by_arm.get(arm)
        if dm is not None and len(dm) > 0:
            ax.scatter(dm[:, 0], dm[:, 1], s=8, alpha=0.45,
                       color=ARM_COLORS.get(arm, "gray"), zorder=3)
        if historical_cloud is not None and n_hist > 0:
            ax.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                       s=8, color="black", zorder=10)
        ax.scatter(*historical_median, marker="*", s=160, color="black", zorder=11)
        _setup_ax(ax, f"({panel_letters[i+1]}) {ARM_LABELS.get(arm, arm)}")

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].set_visible(False)

    # --- Single shared legend below all panels ---
    legend_handles = []
    if n_kirsch > 0:
        legend_handles.append(mlines.Line2D(
            [], [], color="#a0c4e8", marker="o", linestyle="None",
            markersize=5, alpha=0.5, label=f"Kirsch library (n={n_kirsch})"))
    if n_hist > 0:
        legend_handles.append(mlines.Line2D(
            [], [], color="black", marker="o", linestyle="None",
            markersize=5, label=f"Hist. T-blocks (n={n_hist})"))
    legend_handles.append(mlines.Line2D(
        [], [], color="black", marker="*", linestyle="None",
        markersize=9, label="T-block median"))
    for arm in arms:
        dm = pareto_by_arm.get(arm)
        n_arm = len(dm) if dm is not None else 0
        legend_handles.append(mlines.Line2D(
            [], [], color=ARM_COLORS.get(arm, "gray"), marker="o",
            linestyle="None", markersize=6, alpha=0.7,
            label=f"{ARM_LABELS.get(arm, arm)} (n={n_arm})"))

    ncol = min(len(legend_handles), 4)
    fig.suptitle("2D Pareto drought-space by constraint formulation", fontsize=12)
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), ncol=ncol,
               fontsize=8, framealpha=0.9)
    return fig, axes_flat[:n_panels]


def plot_pareto_3d_overlay(
    pareto_by_arm: Dict[str, np.ndarray],
    historical_cloud: Optional[np.ndarray] = None,
    historical_point: Optional[np.ndarray] = None,
    objective_labels: Tuple[str, str, str] = (
        "mean duration", "mean avg severity", "peak severity month"),
    figsize: Tuple[float, float] = (8.0, 6.5),
) -> Tuple[plt.Figure, plt.Axes]:
    """Combined 3D drought-space scatter, one colour per arm (legacy)."""
    apply_style()
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    if historical_cloud is not None and len(historical_cloud) > 0:
        ax.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                   historical_cloud[:, 2],
                   s=10, alpha=0.2, color="gray",
                   label=f"historical (n={len(historical_cloud)})")
    for arm, dm in pareto_by_arm.items():
        if dm is None or len(dm) == 0:
            continue
        ax.scatter(dm[:, 0], dm[:, 1], dm[:, 2],
                   s=14, alpha=0.5,
                   color=ARM_COLORS.get(arm, "gray"),
                   label=f"{ARM_LABELS.get(arm, arm)} (n={len(dm)})")
    if historical_point is not None:
        ax.scatter(historical_point[0], historical_point[1],
                   historical_point[2],
                   marker="*", s=180, color="black",
                   label="historical point")
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_zlabel(objective_labels[2])
    ax.set_title("Pareto drought-space overlay (3D)")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    return fig, ax


def plot_pareto_3d_subpanels(
    pareto_by_arm: Dict[str, np.ndarray],
    historical_cloud: Optional[np.ndarray] = None,
    historical_median: Optional[np.ndarray] = None,
    kirsch_cloud: Optional[np.ndarray] = None,
    objective_labels: Tuple[str, str, str] = (
        "mean duration", "mean avg severity", "peak severity month"),
    arm_order: Optional[List[str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
) -> Tuple[plt.Figure, list]:
    """Multi-panel 3D Pareto, one subplot per arm plus a reference panel.

    Each 3D subplot uses the same axis ranges. Matplotlib 3D axes do not
    support shared-axis links, so limits are set manually. The reference
    panel shows the Kirsch library + historical T-block cloud.
    """
    apply_style()
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    arms = arm_order if arm_order is not None else list(pareto_by_arm.keys())
    n_panels = len(arms) + 1

    ncols = min(3, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    if figsize is None:
        figsize = (5.0 * ncols, 4.5 * nrows)

    fig = plt.figure(figsize=figsize)
    panel_letters = "abcdefghijklmnop"

    # Global limits
    all_x, all_y, all_z = [], [], []
    for src in [historical_cloud, kirsch_cloud]:
        if src is not None and len(src) >= 3:
            all_x.extend(src[:, 0])
            all_y.extend(src[:, 1])
            all_z.extend(src[:, 2])
    for dm in pareto_by_arm.values():
        if dm is not None and len(dm) > 0 and dm.shape[1] >= 3:
            all_x.extend(dm[:, 0])
            all_y.extend(dm[:, 1])
            all_z.extend(dm[:, 2])
    if not all_x:
        return fig, []

    def _lims(vals):
        lo, hi = min(vals), max(vals)
        pad = max(0.05 * (hi - lo), 0.1)
        return lo - pad, hi + pad

    xlim, ylim, zlim = _lims(all_x), _lims(all_y), _lims(all_z)

    def _make_3d(pos, title):
        ax = fig.add_subplot(nrows, ncols, pos, projection="3d")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)
        ax.set_xlabel(objective_labels[0], fontsize=8, labelpad=4)
        ax.set_ylabel(objective_labels[1], fontsize=8, labelpad=4)
        ax.set_zlabel(objective_labels[2], fontsize=8, labelpad=4)
        ax.set_title(title, fontsize=9)
        ax.tick_params(labelsize=7)
        return ax

    n_hist = len(historical_cloud) if historical_cloud is not None else 0
    n_kirsch = len(kirsch_cloud) if kirsch_cloud is not None else 0

    # --- Panel 0: Reference ---
    ax0 = _make_3d(1, f"({panel_letters[0]}) Kirsch + historical T-blocks")
    if kirsch_cloud is not None and n_kirsch > 0 and kirsch_cloud.shape[1] >= 3:
        kc = kirsch_cloud
        if n_kirsch > 1500:
            rng = np.random.default_rng(0)
            kc = kc[rng.choice(n_kirsch, size=1500, replace=False)]
        ax0.scatter(kc[:, 0], kc[:, 1], kc[:, 2],
                    s=5, alpha=0.10, color="#a0c4e8")
    if historical_cloud is not None and n_hist > 0 and historical_cloud.shape[1] >= 3:
        ax0.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                    historical_cloud[:, 2], s=10, color="black", zorder=10)
    if historical_median is not None:
        ax0.scatter(historical_median[0], historical_median[1],
                    historical_median[2], marker="*", s=140, color="black", zorder=11)

    # --- Arm panels ---
    ax_list = [ax0]
    for i, arm in enumerate(arms):
        ax = _make_3d(i + 2, f"({panel_letters[i+1]}) {ARM_LABELS.get(arm, arm)}")
        dm = pareto_by_arm.get(arm)
        if dm is not None and len(dm) > 0 and dm.shape[1] >= 3:
            ax.scatter(dm[:, 0], dm[:, 1], dm[:, 2],
                       s=8, alpha=0.45, color=ARM_COLORS.get(arm, "gray"))
        if historical_cloud is not None and n_hist > 0 and historical_cloud.shape[1] >= 3:
            ax.scatter(historical_cloud[:, 0], historical_cloud[:, 1],
                       historical_cloud[:, 2], s=10, color="black", zorder=10)
        if historical_median is not None:
            ax.scatter(historical_median[0], historical_median[1],
                       historical_median[2], marker="*", s=140,
                       color="black", zorder=11)
        ax_list.append(ax)

    # --- Single shared legend below all panels ---
    legend_handles = []
    if n_kirsch > 0:
        legend_handles.append(mlines.Line2D(
            [], [], color="#a0c4e8", marker="o", linestyle="None",
            markersize=5, alpha=0.5, label=f"Kirsch library (n={n_kirsch})"))
    if n_hist > 0:
        legend_handles.append(mlines.Line2D(
            [], [], color="black", marker="o", linestyle="None",
            markersize=5, label=f"Hist. T-blocks (n={n_hist})"))
    legend_handles.append(mlines.Line2D(
        [], [], color="black", marker="*", linestyle="None",
        markersize=9, label="T-block median"))
    for arm in arms:
        dm = pareto_by_arm.get(arm)
        n_arm = len(dm) if dm is not None else 0
        legend_handles.append(mlines.Line2D(
            [], [], color=ARM_COLORS.get(arm, "gray"), marker="o",
            linestyle="None", markersize=6, alpha=0.7,
            label=f"{ARM_LABELS.get(arm, arm)} (n={n_arm})"))

    ncol = min(len(legend_handles), 4)
    fig.suptitle("3D Pareto drought-space by constraint formulation", fontsize=11)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), ncol=ncol,
               fontsize=8, framealpha=0.9)
    return fig, ax_list


def plot_manhattan_distribution(
    manhattan_by_arm: Dict[str, np.ndarray],
    figsize: Tuple[float, float] = (9, 4),
) -> Tuple[plt.Figure, np.ndarray]:
    """Histogram + empirical CDF of the Manhattan objective per arm.

    Headline falsification figure: if the DV-uniform arm suppresses severe
    droughts, its Manhattan values (distance to anti-ideal) will sit to
    the right of the hydrologic arm's.
    """
    apply_style()
    fig, (ax_h, ax_c) = plt.subplots(1, 2, figsize=figsize)
    xmax = max(
        (np.max(v) for v in manhattan_by_arm.values() if len(v) > 0),
        default=1.0,
    )
    bins = np.linspace(0.0, xmax * 1.05, 40)
    for arm, values in manhattan_by_arm.items():
        values = np.asarray(values, dtype=float)
        if values.size == 0:
            continue
        color = ARM_COLORS.get(arm, "gray")
        label = ARM_LABELS.get(arm, arm)
        med = float(np.median(values))
        ax_h.hist(values, bins=bins, alpha=0.45, color=color,
                  label=f"{label} (n={len(values)}, median={med:.2f})",
                  density=True)
        ax_h.axvline(med, color=color, linestyle="--",
                     linewidth=1.2, alpha=0.9)
        sv = np.sort(values)
        cdf = np.arange(1, len(sv) + 1) / len(sv)
        ax_c.plot(sv, cdf, color=color, linewidth=2, label=label)
    ax_h.set_xlabel(r"Manhattan $f_{K+1}$ (distance to $D^{*}$)")
    ax_h.set_ylabel("density")
    ax_h.set_title("Manhattan-objective distribution")
    ax_h.legend(fontsize=8)
    ax_c.set_xlabel(r"Manhattan $f_{K+1}$")
    ax_c.set_ylabel("empirical CDF")
    ax_c.set_title("Manhattan-objective CDF")
    ax_c.legend(fontsize=8)
    ax_c.set_ylim(0.0, 1.0)
    fig.tight_layout()
    return fig, np.array([ax_h, ax_c])


def plot_dv_distributions(
    dvs_by_arm: Dict[str, np.ndarray],
    figsize: Tuple[float, float] = (9, 4),
    n_dv_sample_rows: Optional[int] = 50,
) -> Tuple[plt.Figure, np.ndarray]:
    """QQ and histogram of flattened Pareto DVs vs ``U[0, 1]`` per arm.

    ``dvs_by_arm[arm]`` is a 2D array ``(n_pareto, n_dvs)``. We visualize
    both the *aggregate* empirical CDF across all Pareto members (left,
    QQ vs uniform) and a marginal histogram (right) so the reader sees
    both the CDF-level deviation (KS-style) and tail clustering.
    """
    apply_style()
    fig, (ax_qq, ax_h) = plt.subplots(1, 2, figsize=figsize)
    ax_qq.plot([0, 1], [0, 1], color="gray", linestyle=":", linewidth=1,
               label="U[0,1] reference")
    for arm, dvs in dvs_by_arm.items():
        if dvs is None or len(dvs) == 0:
            continue
        arr = np.asarray(dvs, dtype=float)
        if arr.ndim == 2 and n_dv_sample_rows and arr.shape[0] > n_dv_sample_rows:
            rng = np.random.default_rng(0)
            pick = rng.choice(arr.shape[0], size=n_dv_sample_rows, replace=False)
            arr = arr[pick]
        flat = arr.reshape(-1)
        color = ARM_COLORS.get(arm, "gray")
        label = ARM_LABELS.get(arm, arm)
        sv = np.sort(flat)
        p = np.arange(1, len(sv) + 1) / (len(sv) + 1)
        ax_qq.plot(p, sv, color=color, linewidth=1.4,
                   label=f"{label} (n={len(flat)})")
        ax_h.hist(flat, bins=40, alpha=0.40, color=color, density=True,
                  label=label)
    ax_qq.set_xlabel("theoretical U[0,1] quantile")
    ax_qq.set_ylabel("empirical DV quantile")
    ax_qq.set_title("Aggregate DV QQ vs U[0,1]")
    ax_qq.legend(fontsize=8)
    ax_qq.set_xlim(0, 1)
    ax_qq.set_ylim(0, 1)
    ax_h.set_xlabel("DV value")
    ax_h.set_ylabel("density")
    ax_h.set_title("DV marginal histogram (aggregate across Pareto)")
    ax_h.axhline(1.0, color="gray", linestyle=":",
                 linewidth=1, label="U[0,1] density")
    ax_h.legend(fontsize=8)
    fig.tight_layout()
    return fig, np.array([ax_qq, ax_h])


def plot_dv_tail_mass(
    dvs_by_arm: Dict[str, np.ndarray],
    tail_bounds: Tuple[float, float] = (0.05, 0.95),
    figsize: Tuple[float, float] = (8.5, 5.0),
) -> Tuple[plt.Figure, np.ndarray]:
    """Compare per-arm fraction of DV mass in the boundary tails.

    For each Pareto member we compute ``frac(dv <= lo) + frac(dv >= hi)``
    — the fraction of DVs landing outside ``tail_bounds``. Under a true
    ``U[0, 1]`` DV vector, this equals ``(lo + 1 - hi)``; for
    ``tail_bounds = (0.05, 0.95)`` the expected value is 0.10. A
    tail-sensitive uniformity constraint should pull the per-trace
    distribution toward that reference line.

    Left panel: histogram of per-trace tail mass, one colour per arm,
    with the uniform reference vertical line. Right panel: empirical
    CDF of the same quantity.
    """
    apply_style()
    fig, (ax_h, ax_c) = plt.subplots(1, 2, figsize=figsize)
    lo, hi = tail_bounds
    reference = lo + (1.0 - hi)
    all_frac = []
    per_arm: Dict[str, np.ndarray] = {}
    for arm, dvs in dvs_by_arm.items():
        if dvs is None or len(dvs) == 0:
            continue
        arr = np.asarray(dvs, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        frac = np.mean((arr <= lo) | (arr >= hi), axis=1)
        per_arm[arm] = frac
        all_frac.append(frac)
    if not per_arm:
        return fig, np.array([ax_h, ax_c])
    xmax = max(float(np.max(v)) for v in per_arm.values())
    bins = np.linspace(0.0, xmax * 1.05 + 1e-6, 40)
    for arm, frac in per_arm.items():
        color = ARM_COLORS.get(arm, "gray")
        label = ARM_LABELS.get(arm, arm)
        ax_h.hist(frac, bins=bins, alpha=0.45, color=color,
                  label=f"{label} (n={len(frac)}, "
                        f"median={float(np.median(frac)):.3f})",
                  density=True)
        ax_h.axvline(float(np.median(frac)), color=color, linestyle="--",
                     linewidth=1.2, alpha=0.9)
        sv = np.sort(frac)
        cdf = np.arange(1, len(sv) + 1) / len(sv)
        ax_c.plot(sv, cdf, color=color, linewidth=2, label=label)
    ax_h.axvline(reference, color="black", linestyle=":",
                 linewidth=1.5,
                 label=f"U[0,1] reference = {reference:.2f}")
    ax_c.axvline(reference, color="black", linestyle=":",
                 linewidth=1.5,
                 label=f"U[0,1] reference = {reference:.2f}")
    ax_h.set_xlabel(f"per-trace fraction of DVs in tails "
                    f"[0, {lo}] ∪ [{hi}, 1]")
    ax_h.set_ylabel("density")
    ax_h.set_title("DV tail-mass distribution per arm")
    ax_h.legend(fontsize=7, loc="best")
    ax_c.set_xlabel("per-trace tail fraction")
    ax_c.set_ylabel("empirical CDF")
    ax_c.set_title("DV tail-mass CDF per arm")
    ax_c.legend(fontsize=8, loc="lower right")
    ax_c.set_ylim(0.0, 1.0)
    fig.tight_layout()
    return fig, np.array([ax_h, ax_c])


def plot_per_trace_stats(
    stats_by_arm: Dict[str, Dict[str, np.ndarray]],
    hist_stats: Dict[str, np.ndarray],
    stat_order: Sequence[str],
    figsize: Tuple[float, float] = (10.5, 5.5),
) -> Tuple[plt.Figure, np.ndarray]:
    """Grid of boxplots — one per summary statistic, one box per arm + historical.

    ``stats_by_arm[arm][stat_name]`` is a 1D array with one value per
    Pareto trace; ``hist_stats[stat_name]`` is the same shape over
    historical blocks. Reviewers can eyeball whether the DV-uniform arm's
    per-trace summaries stay inside the historical envelope (near-
    identical to historical blocks) while the hydrologic arm pushes into
    the drought tail.
    """
    apply_style()
    n = len(stat_order)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = axes.flatten()

    for idx, stat in enumerate(stat_order):
        ax = axes_flat[idx]
        groups = [np.asarray(hist_stats[stat], dtype=float)]
        labels = ["historical"]
        colors = [COLORS["historical"]]
        for arm, arm_stats in stats_by_arm.items():
            vals = arm_stats.get(stat, np.array([]))
            groups.append(np.asarray(vals, dtype=float))
            labels.append(ARM_LABELS.get(arm, arm))
            colors.append(ARM_COLORS.get(arm, "gray"))
        bp = ax.boxplot(groups, labels=labels, patch_artist=True,
                        showfliers=True)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        ax.set_title(stat, fontsize=10)
        ax.tick_params(axis="x", rotation=15, labelsize=8)

    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.tight_layout()
    return fig, axes


def plot_infeasibility_bar(
    infeasibility_by_arm: Dict[str, float],
    figsize: Tuple[float, float] = (5.5, 4.0),
) -> Tuple[plt.Figure, plt.Axes]:
    """Simple bar chart of final infeasible-evaluation fraction per arm."""
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)
    arms = list(infeasibility_by_arm.keys())
    fractions = [infeasibility_by_arm[a] for a in arms]
    colors = [ARM_COLORS.get(a, "gray") for a in arms]
    ax.bar(arms, fractions, color=colors, alpha=0.8)
    for i, f in enumerate(fractions):
        ax.text(i, f, f"{f:.1%}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("infeasible / total evaluations")
    ax.set_title("MOEA infeasibility rate per arm (final)")
    ax.set_ylim(0.0, max(fractions + [0.05]) * 1.25)
    fig.tight_layout()
    return fig, ax
