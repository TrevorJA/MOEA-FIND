"""Drought characteristic space visualization.

Reusable plotting functions for visualizing Pareto fronts and coverage
in drought characteristic space. Supports 2D scatter, 3D projections,
marginal histograms, and coverage comparison panels.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.plotting.style import COLORS, apply_style


def plot_drought_scatter(
    ax: plt.Axes,
    drought_metrics: np.ndarray,
    color: str = COLORS["empirical"],
    label: str = "MOEA-FIND",
    marker_size: float = 8,
    alpha: float = 0.6,
    historical_point: Optional[Tuple[float, float]] = None,
    anti_ideal: Optional[np.ndarray] = None,
    historical_cloud: Optional[np.ndarray] = None,
    historical_cloud_label: str = "Historical T-year blocks",
):
    """Plot drought metrics as scatter on given axes.

    Args:
        ax: Matplotlib axes.
        drought_metrics: Array of shape (n, 2) with the two drought
            characteristic values per Pareto member.
        color: Point color.
        label: Legend label.
        marker_size: Scatter point size.
        alpha: Point transparency.
        historical_point: Optional (duration, intensity) of the
            historical record as a single point (e.g. the
            single-window mean).
        anti_ideal: Optional anti-ideal point to mark.
        historical_cloud: Optional ``(n_blocks, 2)`` array with per-
            historical-block drought characteristics. Each block is
            plotted as a small grey marker so the Pareto archive can
            be compared against the historical *distribution* of
            T-year drought metrics rather than a single point.
        historical_cloud_label: Legend label for the cloud.
    """
    if historical_cloud is not None and len(historical_cloud):
        ax.scatter(
            historical_cloud[:, 0], historical_cloud[:, 1],
            s=18, alpha=0.55, facecolors="none",
            edgecolors="#4a4a4a", linewidths=0.8,
            label=f"{historical_cloud_label} (n={len(historical_cloud)})",
            rasterized=True, zorder=5,
        )

    ax.scatter(drought_metrics[:, 0], drought_metrics[:, 1],
               s=marker_size, alpha=alpha, c=color, label=label, rasterized=True)

    if historical_point is not None:
        ax.plot(historical_point[0], historical_point[1], "k*",
                markersize=12, label="Historical mean", zorder=10)

    if anti_ideal is not None:
        ax.plot(anti_ideal[0], anti_ideal[1], "rx", markersize=10,
                markeredgewidth=2, label="Anti-ideal $D^*$", zorder=10)


def plot_coverage_comparison(
    datasets: List[Tuple[str, np.ndarray, str]],
    historical_point: Optional[Tuple[float, float]] = None,
    anti_ideal: Optional[np.ndarray] = None,
    objective_labels: Tuple[str, str] = ("Mean Duration (months)", "Mean Intensity (cfs)"),
    figsize: Tuple[float, float] = (16, 4),
) -> Tuple[plt.Figure, np.ndarray]:
    """Side-by-side coverage comparison of multiple methods.

    Args:
        datasets: List of (label, metrics_array, color) tuples.
        historical_point: Optional historical mean point.
        anti_ideal: Optional anti-ideal point.
        objective_labels: Axis labels.
        figsize: Figure size.

    Returns:
        (fig, axes) tuple.
    """
    apply_style()
    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=figsize, sharex=True, sharey=True)
    if n == 1:
        axes = [axes]

    for ax, (label, metrics, color) in zip(axes, datasets):
        plot_drought_scatter(ax, metrics, color=color, label=label,
                             historical_point=historical_point,
                             anti_ideal=anti_ideal)
        ax.set_xlabel(objective_labels[0])
        ax.set_title(f"{label} (n={len(metrics)})", fontsize=10)
        ax.legend(fontsize=8, loc="upper left")

    axes[0].set_ylabel(objective_labels[1])
    fig.tight_layout()
    return fig, np.array(axes)


def plot_scatter_with_marginals(
    drought_metrics: np.ndarray,
    color: str = COLORS["empirical"],
    title: str = "",
    historical_point: Optional[Tuple[float, float]] = None,
    anti_ideal: Optional[np.ndarray] = None,
    historical_cloud: Optional[np.ndarray] = None,
    objective_labels: Tuple[str, str] = ("Mean Duration (months)", "Mean Intensity (cfs)"),
    figsize: Tuple[float, float] = (8, 8),
) -> plt.Figure:
    """Scatter plot with marginal histograms.

    Args:
        drought_metrics: Array of shape ``(n, 2)``.
        color: Point/histogram color for the Pareto archive.
        title: Figure title.
        historical_point: Optional historical mean (single-window).
        anti_ideal: Optional anti-ideal.
        historical_cloud: Optional ``(n_blocks, 2)`` per-block
            drought characteristics from
            :func:`src.historical_blocks.compute_historical_block_chars`.
            Drawn as open grey markers in the main panel and as a
            step-histogram on each marginal so the Pareto archive can
            be compared against the historical T-year block
            distribution rather than a single point.
        objective_labels: Axis labels.
        figsize: Figure size.

    Returns:
        Figure object.
    """
    apply_style()
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(4, 4, figure=fig, hspace=0.05, wspace=0.05)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    # Main scatter (Pareto + optional historical cloud + markers)
    plot_drought_scatter(ax_main, drought_metrics, color=color,
                         historical_point=historical_point,
                         anti_ideal=anti_ideal,
                         historical_cloud=historical_cloud)
    ax_main.set_xlabel(objective_labels[0])
    ax_main.set_ylabel(objective_labels[1])
    ax_main.legend(fontsize=8)

    # Top marginal (x axis = objective[0])
    ax_top.hist(drought_metrics[:, 0], bins=30, color=color, alpha=0.6,
                edgecolor="black", linewidth=0.3, label="MOEA-FIND")
    if historical_cloud is not None and len(historical_cloud):
        ax_top.hist(historical_cloud[:, 0], bins=30,
                    histtype="step", color="#4a4a4a", linewidth=1.0,
                    label="Historical blocks")
    if historical_point:
        ax_top.axvline(historical_point[0], color="black", linestyle="--", linewidth=0.8)
    ax_top.tick_params(labelbottom=False)
    ax_top.set_ylabel("Count")

    # Right marginal (y axis = objective[1])
    ax_right.hist(drought_metrics[:, 1], bins=30, orientation="horizontal",
                  color=color, alpha=0.6, edgecolor="black", linewidth=0.3)
    if historical_cloud is not None and len(historical_cloud):
        ax_right.hist(historical_cloud[:, 1], bins=30, orientation="horizontal",
                      histtype="step", color="#4a4a4a", linewidth=1.0)
    if historical_point:
        ax_right.axhline(historical_point[1], color="black", linestyle="--", linewidth=0.8)
    ax_right.tick_params(labelleft=False)
    ax_right.set_xlabel("Count")

    if title:
        fig.suptitle(title, fontsize=12, y=1.01)

    return fig


def plot_drought_space_3d(
    drought_metrics: np.ndarray,
    anti_ideal: np.ndarray,
    objective_labels: Tuple[str, str, str] = ("$D_1$", "$D_2$", "$D_3$"),
    historical_point: Optional[np.ndarray] = None,
    historical_cloud: Optional[np.ndarray] = None,
    title: str = "",
    figsize: Tuple[float, float] = (9, 7),
    include_anti_ideal: bool = False,
):
    """3D scatter of three drought objectives with an optional
    per-block historical cloud.

    The axes are zoomed to the joint range of the Pareto archive,
    historical point, and historical cloud (the anti-ideal is typically
    far outside this range and is annotated in the title rather than
    drawn, so the Pareto cluster remains readable).

    Args:
        drought_metrics: Array of shape ``(n, 3)`` — the three drought
            characteristics Borg saw as objectives, in order.
        anti_ideal: 3-D anti-ideal point. Used for colouring by
            Manhattan distance and reported in the title. Only drawn
            on the axes when ``include_anti_ideal=True``.
        objective_labels: Axis labels in the same order as
            ``drought_metrics`` columns.
        historical_point: Optional 3-D single-window historical point.
        historical_cloud: Optional ``(n_blocks, 3)`` per-block historical
            drought characteristics from
            :func:`src.historical_blocks.compute_historical_block_chars`.
        title: Figure title prefix (appended with hyperplane + D* info).
        figsize: Figure size.
        include_anti_ideal: When True the anti-ideal marker is drawn and
            the axes expanded to contain it — useful once for context,
            but at default headroom the Pareto is typically crushed.

    Returns:
        Matplotlib figure.
    """
    import matplotlib.pyplot as _plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    apply_style()

    dm = np.asarray(drought_metrics, dtype=float)
    D_star = np.asarray(anti_ideal, dtype=float)
    assert dm.shape[1] == 3, "3D scatter requires exactly three objectives"
    assert D_star.shape == (3,), "anti_ideal must be 3-D"

    manh = np.sum(np.abs(dm - D_star), axis=1)

    fig = _plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(dm[:, 0], dm[:, 1], dm[:, 2],
                    c=manh, cmap="viridis_r", s=22, alpha=0.9,
                    edgecolor="none", label="MOEA-FIND")

    if historical_cloud is not None and len(historical_cloud):
        ax.scatter(
            historical_cloud[:, 0], historical_cloud[:, 1], historical_cloud[:, 2],
            facecolors="none", edgecolors="#4a4a4a", s=28, linewidth=0.9,
            label=f"Historical blocks (n={len(historical_cloud)})",
        )

    if historical_point is not None:
        hp = np.asarray(historical_point, dtype=float)
        ax.scatter(*hp, marker="*", color="black", s=280,
                   edgecolor="white", linewidth=0.8,
                   label="Historical mean", zorder=10)

    # Axis ranges: default to Pareto + historical cloud + historical point.
    xs = [dm[:, 0]]
    ys = [dm[:, 1]]
    zs = [dm[:, 2]]
    if historical_cloud is not None and len(historical_cloud):
        xs.append(historical_cloud[:, 0])
        ys.append(historical_cloud[:, 1])
        zs.append(historical_cloud[:, 2])
    if historical_point is not None:
        xs.append(np.array([historical_point[0]]))
        ys.append(np.array([historical_point[1]]))
        zs.append(np.array([historical_point[2]]))
    if include_anti_ideal:
        ax.scatter(*D_star, marker="X", color="red", s=180,
                   edgecolor="black", linewidth=0.6,
                   label=r"anti-ideal $D^*$")
        xs.append(np.array([D_star[0]]))
        ys.append(np.array([D_star[1]]))
        zs.append(np.array([D_star[2]]))
    xs_all = np.concatenate(xs); ys_all = np.concatenate(ys); zs_all = np.concatenate(zs)
    pad_x = 0.05 * (np.ptp(xs_all) + 1e-9)
    pad_y = 0.05 * (np.ptp(ys_all) + 1e-9)
    pad_z = 0.05 * (np.ptp(zs_all) + 1e-9)
    ax.set_xlim(xs_all.min() - pad_x, xs_all.max() + pad_x)
    ax.set_ylim(ys_all.min() - pad_y, ys_all.max() + pad_y)
    ax.set_zlim(zs_all.min() - pad_z, zs_all.max() + pad_z)
    ax.set_xlabel(objective_labels[0])
    ax.set_ylabel(objective_labels[1])
    ax.set_zlabel(objective_labels[2])

    subtitle = (
        f"anti-ideal $D^*=$({D_star[0]:.2f}, {D_star[1]:.2f}, {D_star[2]:.2f})"
        + ("" if include_anti_ideal else " (off-plot)")
    )
    full_title = f"{title}\n{subtitle}" if title else subtitle
    ax.set_title(full_title, fontsize=10)

    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.08)
    cb.set_label(r"Manhattan distance $\|D-D^*\|_1$")
    ax.legend(loc="upper left", fontsize=9)
    return fig


def plot_3d_projections(
    drought_metrics: np.ndarray,
    color: str = COLORS["empirical"],
    objective_labels: Tuple[str, str, str] = ("$D_1$", "$D_2$", "$D_3$"),
    title: str = "3D Pareto Front Projections",
    figsize: Tuple[float, float] = (14, 5),
) -> Tuple[plt.Figure, np.ndarray]:
    """Pairwise 2D projections of 3D drought space.

    Args:
        drought_metrics: Array of shape (n, 3).
        color: Point color.
        objective_labels: Labels for the three objectives.
        title: Figure title.
        figsize: Figure size.

    Returns:
        (fig, axes) tuple.
    """
    apply_style()
    pairs = [(0, 1), (0, 2), (1, 2)]
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    for ax, (i, j) in zip(axes, pairs):
        ax.scatter(drought_metrics[:, i], drought_metrics[:, j],
                   s=4, alpha=0.4, c=color, rasterized=True)
        ax.set_xlabel(objective_labels[i])
        ax.set_ylabel(objective_labels[j])
        ax.set_aspect("equal")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig, axes


# =============================================================================
# Manuscript Figure 6 — three-way drought-space composite
# =============================================================================
#
# Building blocks used to assemble a single Figure 6 that replaces the three
# redundant drought-space figures produced by exp04 / exp11. The composite
# shows a 2D scatter + marginals, a 3D scatter, and four sample monthly
# flow traces (2 historical blocks + 2 nearest-neighbor MOEA-FIND Pareto
# traces) with SSI-based drought events shaded.


def plot_three_way_scatter_marginals(
    moea: np.ndarray,
    kirsch: np.ndarray,
    hist_blocks: np.ndarray,
    ax_main: plt.Axes,
    ax_top: plt.Axes,
    ax_right: plt.Axes,
    anti_ideal: Optional[Tuple[float, float]] = None,
    historical_point: Optional[Tuple[float, float]] = None,
    objective_labels: Tuple[str, str] = (
        "mean duration (months)", "mean avg severity (SSI)",
    ),
    bins: int = 32,
    highlights: Optional[List[Dict]] = None,
) -> None:
    """2D scatter + marginal step histograms for three populations.

    Main axes show MOEA-FIND Pareto (red filled), Kirsch random ensemble
    (gray filled, low alpha, drawn underneath), and historical T-year
    blocks (black open circles). Marginals are three stacked step
    histograms (solid/dashed/dotted) so overlapping distributions stay
    legible.
    """
    moea = np.asarray(moea, dtype=float)
    kirsch = np.asarray(kirsch, dtype=float)
    hist_blocks = np.asarray(hist_blocks, dtype=float)

    c_moea = COLORS["parametric"]
    c_kirsch = COLORS["muted"]
    c_hist = COLORS["historical"]

    # -- Main scatter: Kirsch first (background), MOEA next, historical on top --
    if len(kirsch):
        ax_main.scatter(
            kirsch[:, 0], kirsch[:, 1],
            s=2.5, alpha=0.12, c=c_kirsch, edgecolor="none",
            rasterized=True, zorder=1,
            label=f"Kirsch random (n={len(kirsch)})",
        )
    if len(moea):
        ax_main.scatter(
            moea[:, 0], moea[:, 1],
            s=4, alpha=0.35, c=c_moea, edgecolor="none",
            rasterized=True, zorder=3,
            label=f"MOEA-FIND (n={len(moea)})",
        )
    if len(hist_blocks):
        ax_main.scatter(
            hist_blocks[:, 0], hist_blocks[:, 1],
            s=26, facecolors="none", edgecolors=c_hist, linewidths=1.0,
            zorder=5,
            label=f"historical blocks (n={len(hist_blocks)})",
        )
    if historical_point is not None:
        ax_main.plot(historical_point[0], historical_point[1], "*",
                     color=c_hist, markersize=14,
                     markeredgecolor="white", markeredgewidth=0.8,
                     label="historical mean", zorder=10)

    # Highlight the focal-scenario points used by the timeseries panels.
    # Each highlight dict: {point: (x,y), color, marker, label}.
    if highlights:
        for h in highlights:
            x, y = h["point"]
            mk = h.get("marker", "o")
            col = h.get("color", c_hist)
            lbl = h.get("label", "")
            ax_main.plot(
                x, y, marker=mk, markersize=11,
                markerfacecolor="none", markeredgecolor=col,
                markeredgewidth=1.8, linestyle="none", zorder=11,
            )
            ax_main.annotate(
                lbl, xy=(x, y), xytext=(6, 6),
                textcoords="offset points", fontsize=8, fontweight="bold",
                color=col, zorder=12,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          alpha=0.85, edgecolor="none"),
            )

    # Zoom to the data cluster; keep D* off-plot with a corner annotation
    # so the scatter stays readable (anti-ideal often sits far outside).
    xs = np.concatenate([
        a[:, 0] for a in (moea, kirsch, hist_blocks) if len(a)
    ])
    ys = np.concatenate([
        a[:, 1] for a in (moea, kirsch, hist_blocks) if len(a)
    ])
    pad_x = 0.05 * max(np.ptp(xs), 1e-6)
    pad_y = 0.05 * max(np.ptp(ys), 1e-6)
    ax_main.set_xlim(xs.min() - pad_x, xs.max() + pad_x)
    ax_main.set_ylim(ys.min() - pad_y, ys.max() + pad_y)

    if anti_ideal is not None:
        ax_main.text(
            0.98, 0.98,
            rf"$D^\star = ({anti_ideal[0]:.1f},\,{anti_ideal[1]:.1f})$"
            "\n(off-plot)",
            transform=ax_main.transAxes, ha="right", va="top",
            fontsize=7, color=COLORS["anti_ideal"],
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      alpha=0.9, edgecolor=COLORS["anti_ideal"],
                      linewidth=0.6),
        )

    ax_main.set_xlabel(objective_labels[0])
    ax_main.set_ylabel(objective_labels[1])
    ax_main.legend(fontsize=7, loc="upper left", framealpha=0.9)

    # Shared bin edges so the three histograms on each margin align exactly.
    stack_x = np.concatenate([
        a[:, 0] for a in (kirsch, moea, hist_blocks) if len(a)
    ])
    stack_y = np.concatenate([
        a[:, 1] for a in (kirsch, moea, hist_blocks) if len(a)
    ])
    edges_x = np.linspace(stack_x.min(), stack_x.max(), bins + 1)
    edges_y = np.linspace(stack_y.min(), stack_y.max(), bins + 1)

    def _step(ax, values, edges, color, linestyle, label, orientation):
        if not len(values):
            return
        ax.hist(values, bins=edges, histtype="step", density=True,
                color=color, linestyle=linestyle, linewidth=1.2,
                orientation=orientation, label=label)

    # Top marginal (axis 0)
    _step(ax_top, kirsch[:, 0], edges_x, c_kirsch, ":", "Kirsch", "vertical")
    _step(ax_top, moea[:, 0], edges_x, c_moea, "-", "MOEA-FIND", "vertical")
    _step(ax_top, hist_blocks[:, 0], edges_x, c_hist, "--", "historical",
          "vertical")
    ax_top.tick_params(labelbottom=False, left=False, labelleft=False)
    ax_top.spines["left"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["top"].set_visible(False)

    # Right marginal (axis 1)
    _step(ax_right, kirsch[:, 1], edges_y, c_kirsch, ":", None, "horizontal")
    _step(ax_right, moea[:, 1], edges_y, c_moea, "-", None, "horizontal")
    _step(ax_right, hist_blocks[:, 1], edges_y, c_hist, "--", None,
          "horizontal")
    ax_right.tick_params(labelleft=False, bottom=False, labelbottom=False)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["bottom"].set_visible(False)


def plot_three_way_drought_space_3d(
    ax,
    moea: np.ndarray,
    kirsch: np.ndarray,
    hist_blocks: np.ndarray,
    objective_labels: Tuple[str, str, str] = ("$D_1$", "$D_2$", "$D_3$"),
    historical_point: Optional[np.ndarray] = None,
    anti_ideal: Optional[np.ndarray] = None,
    view: Tuple[float, float] = (22.0, -60.0),
) -> None:
    """3D scatter on a provided Axes3D for three populations.

    Kirsch is drawn as a low-alpha gray cloud beneath MOEA-FIND (red).
    Historical blocks are open black circles on top.
    """
    moea = np.asarray(moea, dtype=float)
    kirsch = np.asarray(kirsch, dtype=float)
    hist_blocks = np.asarray(hist_blocks, dtype=float)

    c_moea = COLORS["parametric"]
    c_kirsch = COLORS["muted"]
    c_hist = COLORS["historical"]

    if len(kirsch):
        ax.scatter(kirsch[:, 0], kirsch[:, 1], kirsch[:, 2],
                   s=2.5, alpha=0.06, c=c_kirsch, edgecolor="none",
                   label=f"Kirsch (n={len(kirsch)})")
    if len(moea):
        ax.scatter(moea[:, 0], moea[:, 1], moea[:, 2],
                   s=5, alpha=0.25, c=c_moea, edgecolor="none",
                   label=f"MOEA-FIND (n={len(moea)})")
    if len(hist_blocks):
        ax.scatter(hist_blocks[:, 0], hist_blocks[:, 1], hist_blocks[:, 2],
                   s=24, facecolors="none", edgecolors=c_hist, linewidths=1.0,
                   label=f"historical (n={len(hist_blocks)})")
    if historical_point is not None:
        hp = np.asarray(historical_point, dtype=float)
        ax.scatter(*hp, marker="*", color=c_hist, s=200,
                   edgecolor="white", linewidth=0.8,
                   label="historical mean", zorder=10)

    # Zoom to the data cluster; anti-ideal is reported in the axes title
    # (drawn by the caller) rather than the plot itself because D* is
    # usually far outside the cluster and would crush everything.
    xs = np.concatenate([
        a[:, 0] for a in (moea, kirsch, hist_blocks) if len(a)
    ])
    ys = np.concatenate([
        a[:, 1] for a in (moea, kirsch, hist_blocks) if len(a)
    ])
    zs = np.concatenate([
        a[:, 2] for a in (moea, kirsch, hist_blocks) if len(a)
    ])
    pad_x = 0.05 * max(np.ptp(xs), 1e-6)
    pad_y = 0.05 * max(np.ptp(ys), 1e-6)
    pad_z = 0.05 * max(np.ptp(zs), 1e-6)
    ax.set_xlim(xs.min() - pad_x, xs.max() + pad_x)
    ax.set_ylim(ys.min() - pad_y, ys.max() + pad_y)
    ax.set_zlim(zs.min() - pad_z, zs.max() + pad_z)

    ax.set_xlabel(objective_labels[0], labelpad=-2)
    ax.set_ylabel(objective_labels[1], labelpad=-2)
    ax.set_zlabel(objective_labels[2], labelpad=-2)
    ax.view_init(elev=view[0], azim=view[1])
    for lbl in (ax.get_xticklabels() + ax.get_yticklabels()
                + ax.get_zticklabels()):
        lbl.set_fontsize(7)
    if anti_ideal is not None:
        D = np.asarray(anti_ideal, dtype=float)
        ax.text2D(
            0.02, 0.02,
            rf"$D^\star=({D[0]:.1f},\,{D[1]:.1f},\,{D[2]:.1f})$ (off-plot)",
            transform=ax.transAxes, fontsize=7,
            color=COLORS["anti_ideal"],
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      alpha=0.9, edgecolor=COLORS["anti_ideal"],
                      linewidth=0.6),
        )
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)


def _select_extremes(
    hist_block_chars: np.ndarray,
    pareto_metrics: np.ndarray,
    dur_idx: int = 0,
    sev_idx: int = 1,
    percentile: float = 90.0,
) -> Dict[str, int]:
    """Pick representative indices for the Fig 6 timeseries panels.

    Instead of picking the single hard-max on each axis (which can be a
    lone outlier), this picks the member **closest to the P-th percentile**
    on each axis (default P=90). The result is a "strong-but-not-extreme"
    representative that reads better in the manuscript figure.

    Selection rule:
      - Hist-A = historical block closest to the P-th percentile of severity.
      - Hist-B = historical block closest to the P-th percentile of duration.
      - MOEA-A = Pareto member closest to the P-th percentile of severity.
      - MOEA-B = Pareto member closest to the P-th percentile of duration.

    If Hist-A / Hist-B (or MOEA-A / MOEA-B) collide, the second pick
    falls back to the next-closest candidate so every panel shows a
    distinct hydrograph.

    Returns a dict with keys ``{"hist_A", "hist_B", "moea_A", "moea_B"}``.
    """
    h = np.asarray(hist_block_chars, dtype=float)
    p = np.asarray(pareto_metrics, dtype=float)
    assert h.shape[1] >= 2 and p.shape[1] >= 2

    def _pick_near_pct(values: np.ndarray, exclude: Optional[int] = None) -> int:
        target = float(np.percentile(values, percentile))
        order = np.argsort(np.abs(values - target))
        for cand in order:
            if exclude is None or int(cand) != exclude:
                return int(cand)
        return int(order[0])

    i_hist_A = _pick_near_pct(h[:, sev_idx])
    i_hist_B = _pick_near_pct(h[:, dur_idx], exclude=i_hist_A)
    i_moea_A = _pick_near_pct(p[:, sev_idx])
    i_moea_B = _pick_near_pct(p[:, dur_idx], exclude=i_moea_A)

    return {
        "hist_A": i_hist_A,
        "hist_B": i_hist_B,
        "moea_A": i_moea_A,
        "moea_B": i_moea_B,
    }


def plot_trace_with_ssi_events(
    ax: plt.Axes,
    monthly_flow: np.ndarray,
    ssi_calc,
    start_date: Optional[str] = None,
    title: str = "",
    color: str = None,
    threshold_pct: float = 20.0,
):
    """Render a monthly hydrograph with SSI-based drought events shaded.

    Pipeline:
      1. Wrap ``monthly_flow`` as a pandas Series with a datetime index
         starting at ``start_date`` and transform it to SSI with the
         supplied (pre-fitted) SSI calculator — so every panel uses the
         same reference gamma fit.
      2. Detect critical drought events via
         :func:`src.plotting.trace_diagnostics._detect_drought_events`.
      3. Draw the monthly flow line, a horizontal dashed threshold at the
         trace's Pth percentile (default P20) for orientation, and a
         shaded ``axvspan`` over each detected critical event.

    If ``start_date`` is None, the x-axis is rendered as a simple
    water-year index (year 1, 2, ..., T) rather than calendar dates —
    appropriate for synthetic traces with no meaningful calendar year.

    Returns the list of event dicts.
    """
    import pandas as pd
    from src.plotting.trace_diagnostics import _detect_drought_events
    from src.objectives import flows_to_series

    if color is None:
        color = COLORS["empirical"]

    flow = np.asarray(monthly_flow, dtype=float).ravel()
    # SSI calibration depends only on calendar-month seasonality, so we
    # always attach a datetime index for the transform; the *display*
    # axis is separate.
    calibration_start = start_date if start_date else "2100-10-01"
    series = flows_to_series(flow, start_date=calibration_start)
    ssi = ssi_calc.transform(series)
    if hasattr(ssi, "columns"):
        ssi = ssi.iloc[:, 0]
    events = _detect_drought_events(pd.Series(np.asarray(ssi), index=ssi.index))

    # x-axis: real dates if a caller-supplied start, else water-year index.
    if start_date is not None:
        t = series.index
        to_x = lambda ts: ts
    else:
        # Convert datetime → water-year-relative month index (1..N)
        n = len(flow)
        t = np.arange(1, n + 1) / 12.0  # years
        series_dt = series.index
        def to_x(ts):
            return float((ts - series_dt[0]).days) / 365.25 + 1.0 / 12.0

    # Log-y flow is standard practice for drought/low-flow viz; clip at
    # 1 cfs so zero / near-zero months don't blow out the bottom.
    flow_plot = np.clip(flow, 1.0, None)
    ax.plot(t, flow_plot, color=color, linewidth=0.8, alpha=0.9)
    thresh = float(np.percentile(flow, threshold_pct))
    ax.axhline(max(thresh, 1.0), color=COLORS["muted"], linestyle="--",
               linewidth=0.6, alpha=0.7,
               label=f"P{int(threshold_pct)} threshold")
    for ev in events:
        ax.axvspan(to_x(ev["start"]), to_x(ev["end"]), alpha=0.25,
                   color=COLORS["highlight"], zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(bottom=1.0)
    ax.set_ylabel("flow (cfs)")
    if start_date is None:
        ax.set_xlabel("water year in block")
    if title:
        ax.set_title(title, fontsize=9)
    ax.tick_params(axis="x", labelsize=7)
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(25)
        lbl.set_ha("right")
    ax.margins(x=0.01)
    n_ev = len(events)
    ax.text(0.02, 0.95, f"{n_ev} SSI event{'s' if n_ev != 1 else ''}",
            transform=ax.transAxes, fontsize=7, va="top",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      alpha=0.85, edgecolor="none"))
    return events


def plot_fig6_composite(
    pareto_metrics: np.ndarray,
    pareto_traces_1d: List[np.ndarray],
    kirsch_2d: np.ndarray,
    kirsch_3d: Optional[np.ndarray],
    hist_block_chars_2d: np.ndarray,
    hist_block_chars_3d: Optional[np.ndarray],
    hist_blocks_1d: List[np.ndarray],
    ssi_calc,
    anti_ideal: np.ndarray,
    historical_point_2d: Optional[Tuple[float, float]] = None,
    historical_point_3d: Optional[np.ndarray] = None,
    objective_labels_2d: Tuple[str, str] = (
        "mean duration (months)", "mean avg severity (SSI)",
    ),
    objective_labels_3d: Tuple[str, str, str] = (
        "duration", "avg severity", "peak month",
    ),
    historical_start_years: Optional[List[int]] = None,
    figsize: Tuple[float, float] = (13.0, 11.5),
) -> plt.Figure:
    """Assemble manuscript Figure 6.

    Layout:
      Top: (a) 2D scatter + marginals (left) and (b) 3D scatter (right),
      both showing three populations — MOEA-FIND Pareto, Kirsch random
      ensemble, historical T-year blocks.

      Below, stacked vertically spanning the full figure width:
        (c) historical block with max severity
        (d) historical block with max duration
        (e) MOEA-FIND Pareto member with max severity
        (f) MOEA-FIND Pareto member with max duration

    Each timeseries panel shades SSI critical-drought event spans using
    the pre-fitted historical SSI for a consistent reference.
    """
    apply_style()

    pm = np.asarray(pareto_metrics, dtype=float)
    km = np.asarray(kirsch_2d, dtype=float)
    hm = np.asarray(hist_block_chars_2d, dtype=float)
    assert pm.shape[1] >= 2

    picks = _select_extremes(hm, pm, dur_idx=0, sev_idx=1)

    fig = plt.figure(figsize=figsize)
    # Scatters row (~3 units tall) + 4 stacked timeseries rows (~1 unit each).
    outer = GridSpec(
        5, 2, figure=fig,
        height_ratios=[3.0, 1.0, 1.0, 1.0, 1.0],
        hspace=0.65, wspace=0.22,
    )

    # -- (a) 2D scatter with marginals --
    a_gs = outer[0, 0].subgridspec(4, 4, hspace=0.05, wspace=0.05)
    ax_top_a = fig.add_subplot(a_gs[0, 0:3])
    ax_main_a = fig.add_subplot(a_gs[1:4, 0:3], sharex=ax_top_a)
    ax_right_a = fig.add_subplot(a_gs[1:4, 3], sharey=ax_main_a)

    # Build the focal-scenario markers tied to the four timeseries panels.
    c_hist = COLORS["historical"]
    c_moea = COLORS["parametric"]
    highlights = [
        {"point": (float(hm[picks["hist_A"], 0]),
                   float(hm[picks["hist_A"], 1])),
         "color": c_hist, "marker": "s", "label": "(c)"},
        {"point": (float(hm[picks["hist_B"], 0]),
                   float(hm[picks["hist_B"], 1])),
         "color": c_hist, "marker": "D", "label": "(d)"},
        {"point": (float(pm[picks["moea_A"], 0]),
                   float(pm[picks["moea_A"], 1])),
         "color": c_moea, "marker": "s", "label": "(e)"},
        {"point": (float(pm[picks["moea_B"], 0]),
                   float(pm[picks["moea_B"], 1])),
         "color": c_moea, "marker": "D", "label": "(f)"},
    ]

    plot_three_way_scatter_marginals(
        moea=pm[:, :2],
        kirsch=km[:, :2],
        hist_blocks=hm[:, :2],
        ax_main=ax_main_a, ax_top=ax_top_a, ax_right=ax_right_a,
        anti_ideal=(
            (float(anti_ideal[0]), float(anti_ideal[1]))
            if anti_ideal is not None else None
        ),
        historical_point=historical_point_2d,
        objective_labels=objective_labels_2d,
        highlights=highlights,
    )
    ax_top_a.set_title("(a) 2D drought space", fontsize=10, loc="left")

    # -- (b) 3D scatter --
    ax_b = fig.add_subplot(outer[0, 1], projection="3d")
    can_3d = (
        pm.shape[1] >= 3
        and hist_block_chars_3d is not None
        and np.asarray(hist_block_chars_3d).shape[1] >= 3
    )
    if can_3d:
        h3 = np.asarray(hist_block_chars_3d, dtype=float)[:, :3]
        p3 = pm[:, :3]
        k3 = (
            np.asarray(kirsch_3d, dtype=float)[:, :3]
            if kirsch_3d is not None and len(kirsch_3d) else np.empty((0, 3))
        )
        plot_three_way_drought_space_3d(
            ax_b,
            moea=p3, kirsch=k3, hist_blocks=h3,
            objective_labels=objective_labels_3d,
            historical_point=(
                np.asarray(historical_point_3d, dtype=float)
                if historical_point_3d is not None else None
            ),
            anti_ideal=(
                np.asarray(anti_ideal, dtype=float)[:3]
                if anti_ideal is not None else None
            ),
        )
    else:
        ax_b.text(0.5, 0.5, "(3D data unavailable)",
                  transform=ax_b.transAxes, ha="center", va="center")
    ax_b.set_title("(b) 3D drought space", fontsize=10)

    # -- Stacked timeseries panels --
    def _safe_start(i_hist: int) -> Optional[str]:
        if historical_start_years and i_hist < len(historical_start_years):
            return f"{int(historical_start_years[i_hist])}-10-01"
        return None

    def _year_range(i_hist: int) -> str:
        if historical_start_years and i_hist < len(historical_start_years):
            y0 = int(historical_start_years[i_hist])
            y1 = y0 + len(hist_blocks_1d[0]) // 12 - 1
            return f" (WY {y0}–{y1})"
        return ""

    def _chars_annot(row: np.ndarray) -> str:
        parts = []
        if row.size > 0:
            parts.append(f"dur={row[0]:.2f} mo")
        if row.size > 1:
            parts.append(f"sev={row[1]:.2f}")
        if row.size > 2:
            parts.append(f"peak mo={row[2]:.1f}")
        return "  |  ".join(parts)

    i_hist_A = picks["hist_A"]
    i_hist_B = picks["hist_B"]
    i_moea_A = picks["moea_A"]
    i_moea_B = picks["moea_B"]

    panels = [
        (outer[1, :], hist_blocks_1d[i_hist_A],
         f"(c) historical, max-severity block{_year_range(i_hist_A)}  —  "
         f"{_chars_annot(hm[i_hist_A])}",
         COLORS["historical"], _safe_start(i_hist_A)),
        (outer[2, :], hist_blocks_1d[i_hist_B],
         f"(d) historical, max-duration block{_year_range(i_hist_B)}  —  "
         f"{_chars_annot(hm[i_hist_B])}",
         COLORS["historical"], _safe_start(i_hist_B)),
        (outer[3, :], np.asarray(pareto_traces_1d[i_moea_A]),
         f"(e) MOEA-FIND max-severity Pareto member #{i_moea_A}  —  "
         f"{_chars_annot(pm[i_moea_A])}",
         COLORS["parametric"], None),
        (outer[4, :], np.asarray(pareto_traces_1d[i_moea_B]),
         f"(f) MOEA-FIND max-duration Pareto member #{i_moea_B}  —  "
         f"{_chars_annot(pm[i_moea_B])}",
         COLORS["parametric"], None),
    ]
    last_idx = len(panels) - 1
    for i, (spec, flow, title, color, start) in enumerate(panels):
        ax = fig.add_subplot(spec)
        plot_trace_with_ssi_events(
            ax, flow, ssi_calc,
            start_date=start, title=title, color=color,
        )
        # Only show the x-axis label on the bottom-most panel so the
        # stacked layout doesn't collide labels with the next title.
        if i != last_idx:
            ax.set_xlabel("")

    return fig


def plot_epsilon_box_filling(
    drought_metrics: np.ndarray,
    epsilons: np.ndarray,
    lb: np.ndarray,
    ub: np.ndarray,
    figsize: Tuple[float, float] = (8, 6),
) -> Tuple[plt.Figure, plt.Axes]:
    """Visualize epsilon-box grid filling in 2D drought space.

    Shows which epsilon-boxes are occupied by Pareto solutions and
    which are empty. Reports box filling rate.

    Args:
        drought_metrics: Array of shape (n, 2).
        epsilons: Epsilon values for each objective (2,).
        lb: Lower bounds of drought space.
        ub: Upper bounds of drought space.
        figsize: Figure size.

    Returns:
        (fig, ax) tuple.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Compute grid
    n_boxes_x = int(np.ceil((ub[0] - lb[0]) / epsilons[0]))
    n_boxes_y = int(np.ceil((ub[1] - lb[1]) / epsilons[1]))
    total_boxes = n_boxes_x * n_boxes_y

    # Assign points to boxes
    box_x = np.clip(((drought_metrics[:, 0] - lb[0]) / epsilons[0]).astype(int),
                     0, n_boxes_x - 1)
    box_y = np.clip(((drought_metrics[:, 1] - lb[1]) / epsilons[1]).astype(int),
                     0, n_boxes_y - 1)
    occupied = set(zip(box_x, box_y))
    n_occupied = len(occupied)
    fill_rate = n_occupied / total_boxes

    # Draw grid
    for i in range(n_boxes_x + 1):
        x = lb[0] + i * epsilons[0]
        ax.axvline(x, color="lightgray", linewidth=0.3)
    for j in range(n_boxes_y + 1):
        y = lb[1] + j * epsilons[1]
        ax.axhline(y, color="lightgray", linewidth=0.3)

    # Color occupied boxes
    for bx, by in occupied:
        rect = plt.Rectangle(
            (lb[0] + bx * epsilons[0], lb[1] + by * epsilons[1]),
            epsilons[0], epsilons[1],
            facecolor=COLORS["empirical"], alpha=0.15,
        )
        ax.add_patch(rect)

    # Scatter points
    ax.scatter(drought_metrics[:, 0], drought_metrics[:, 1],
               s=6, alpha=0.7, c=COLORS["empirical"], zorder=5)

    ax.set_xlim(lb[0], ub[0])
    ax.set_ylim(lb[1], ub[1])
    ax.set_title(f"Epsilon-Box Filling: {n_occupied}/{total_boxes} "
                 f"({fill_rate:.1%})", fontsize=11)
    ax.set_xlabel("Mean Duration (months)")
    ax.set_ylabel("Mean Intensity (cfs)")

    fig.tight_layout()
    return fig, ax
