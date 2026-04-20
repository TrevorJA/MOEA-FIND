"""Plotting composition helpers for the DV-uniformity ablation SI figures.

Wraps the existing primitives in :mod:`src.plotting.trace_diagnostics` so the
compare script (``workflows/experiments/14_dv_uniformity_compare.py``) can
render side-by-side or overlaid comparisons of the hydrologic arm and the
DV-uniform arm against the same historical block envelope, without
duplicating the envelope-percentile machinery.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from src.plotting.style import COLORS, WATER_YEAR_MONTHS, apply_style
from src.plotting.trace_diagnostics import _acf, _envelope, _fdc


ARM_COLORS: Dict[str, str] = {
    "hydrologic": "#2b6cb0",  # steel blue
    "dv_uniform": "#c05621",  # terracotta
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
                            label=f"{arm} (n={len(traces)}, 10-90%)")
        ax_acf.plot(lags, s_mid, "--", color=ARM_COLORS.get(arm, "gray"),
                    linewidth=1.5, label=f"{arm} median")
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
                            label=f"{arm} (n={len(traces)}, 10-90%)")
        ax_fdc.semilogy(grid_pct, s_mid, "--",
                        color=ARM_COLORS.get(arm, "gray"),
                        linewidth=1.5, label=f"{arm} median")
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
                    linewidth=1.5, label=f"{arm} median")
        ax.set_xticks(months)
        ax.set_xticklabels(WATER_YEAR_MONTHS, rotation=45, fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best", framealpha=0.9)

    fig.tight_layout()
    return fig, axes


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
        ax_h.hist(values, bins=bins, alpha=0.55, color=color,
                  label=f"{arm} (n={len(values)}, "
                        f"median={np.median(values):.2f})",
                  density=True)
        sv = np.sort(values)
        cdf = np.arange(1, len(sv) + 1) / len(sv)
        ax_c.plot(sv, cdf, color=color, linewidth=2, label=arm)
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
        sv = np.sort(flat)
        p = np.arange(1, len(sv) + 1) / (len(sv) + 1)
        ax_qq.plot(p, sv, color=color, linewidth=1.2,
                   label=f"{arm} (n={len(flat)})")
        ax_h.hist(flat, bins=40, alpha=0.45, color=color, density=True,
                  label=f"{arm}")
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
        for arm in ("hydrologic", "dv_uniform"):
            vals = stats_by_arm.get(arm, {}).get(stat, np.array([]))
            groups.append(np.asarray(vals, dtype=float))
            labels.append(arm)
            colors.append(ARM_COLORS[arm])
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
