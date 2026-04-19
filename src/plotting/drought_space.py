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
