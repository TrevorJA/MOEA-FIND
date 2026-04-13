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
):
    """Plot drought metrics as scatter on given axes.

    Args:
        ax: Matplotlib axes.
        drought_metrics: Array of shape (n, 2) with (duration, intensity).
        color: Point color.
        label: Legend label.
        marker_size: Scatter point size.
        alpha: Point transparency.
        historical_point: Optional (duration, intensity) of historical mean.
        anti_ideal: Optional anti-ideal point to mark.
    """
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
    objective_labels: Tuple[str, str] = ("Mean Duration (months)", "Mean Intensity (cfs)"),
    figsize: Tuple[float, float] = (8, 8),
) -> plt.Figure:
    """Scatter plot with marginal histograms.

    Args:
        drought_metrics: Array of shape (n, 2).
        color: Point/histogram color.
        title: Figure title.
        historical_point: Optional historical mean.
        anti_ideal: Optional anti-ideal.
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

    # Main scatter
    plot_drought_scatter(ax_main, drought_metrics, color=color,
                         historical_point=historical_point,
                         anti_ideal=anti_ideal)
    ax_main.set_xlabel(objective_labels[0])
    ax_main.set_ylabel(objective_labels[1])
    ax_main.legend(fontsize=8)

    # Top marginal
    ax_top.hist(drought_metrics[:, 0], bins=30, color=color, alpha=0.6,
                edgecolor="black", linewidth=0.3)
    if historical_point:
        ax_top.axvline(historical_point[0], color="black", linestyle="--", linewidth=0.8)
    ax_top.tick_params(labelbottom=False)
    ax_top.set_ylabel("Count")

    # Right marginal
    ax_right.hist(drought_metrics[:, 1], bins=30, orientation="horizontal",
                  color=color, alpha=0.6, edgecolor="black", linewidth=0.3)
    if historical_point:
        ax_right.axhline(historical_point[1], color="black", linestyle="--", linewidth=0.8)
    ax_right.tick_params(labelleft=False)
    ax_right.set_xlabel("Count")

    if title:
        fig.suptitle(title, fontsize=12, y=1.01)

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
