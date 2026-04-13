"""Coverage-comparison figures (manuscript §6.3, Fig 7).

The headline figure of the paper: MOEA-FIND versus library-and-subsample
baselines in drought-outcome space. Functions here consume the outputs of
scripts 04 (MOEA-FIND Pareto), 05 (Kirsch library), and 06 (library
subsamples) and produce a single composite panel used as Figure 7 plus
the companion SI table.

Sizing conventions: WRR double-column = 7.0 in wide.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.analysis import coverage_metrics
from src.plotting.style import COLORS, apply_style

FeatureSet = Tuple[str, np.ndarray, str]  # (label, points (n,d), color)


def _default_colors_for_labels(labels: Sequence[str]) -> List[str]:
    palette = {
        "MOEA-FIND": COLORS["empirical"],
        "Library-LHS": COLORS["lhs"],
        "Library-Sobol": COLORS["sobol"],
        "Library-Random": COLORS["random"],
        "Full library": COLORS["muted"],
    }
    return [palette.get(l, COLORS["highlight"]) for l in labels]


def fig7_coverage_comparison(
    feature_sets: List[FeatureSet],
    historical_point: Optional[Tuple[float, float]] = None,
    anti_ideal: Optional[np.ndarray] = None,
    feature_labels: Tuple[str, str] = ("Mean duration (months)",
                                       "Mean avg. severity (SSI units)"),
    figsize: Tuple[float, float] = (7.0, 3.6),
) -> plt.Figure:
    """Manuscript Figure 7 — coverage comparison across methods.

    Panels:
        (a) Scatter overlay of MOEA-FIND Pareto vs library subsamples.
        (b) L2-star / NN_CV bar chart (same method ordering).

    Args:
        feature_sets: [(label, (n, 2) array, color)]. The first entry is
            treated as MOEA-FIND and plotted on top.
        historical_point: Drought characteristics of the historical record.
        anti_ideal: Anti-ideal D* used by the Manhattan-norm formulation.
        feature_labels: Axis labels for the drought characteristic space.
        figsize: Figure size in inches.
    """
    apply_style()
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.4, 1.0], wspace=0.3)

    # -- (a) Scatter overlay --
    ax_a = fig.add_subplot(gs[0, 0])
    for label, pts, color in reversed(feature_sets):
        ax_a.scatter(pts[:, 0], pts[:, 1], s=8, alpha=0.55, c=color,
                     label=f"{label} (n={len(pts)})", rasterized=True)
    if historical_point is not None:
        ax_a.plot(historical_point[0], historical_point[1], "k*", ms=10,
                  label="Historical")
    if anti_ideal is not None:
        ax_a.plot(anti_ideal[0], anti_ideal[1], "x",
                  color=COLORS["anti_ideal"], ms=9, mew=2,
                  label=r"Anti-ideal $D^*$")
    ax_a.set_xlabel(feature_labels[0])
    ax_a.set_ylabel(feature_labels[1])
    ax_a.set_title("(a) Coverage in drought space", fontsize=10)
    ax_a.legend(fontsize=7, loc="best", framealpha=0.9)

    # -- (b) Coverage-metric bar chart --
    ax_b = fig.add_subplot(gs[0, 1])
    labels = [fs[0] for fs in feature_sets]
    lb = np.min(np.vstack([fs[1] for fs in feature_sets]), axis=0)
    ub = np.max(np.vstack([fs[1] for fs in feature_sets]), axis=0)
    l2 = []
    nn = []
    for _, pts, _ in feature_sets:
        m = coverage_metrics(pts, lb, ub)
        l2.append(m["L2_star_discrepancy"])
        nn.append(m["nn_cv"])

    x = np.arange(len(labels))
    w = 0.38
    colors = [fs[2] for fs in feature_sets]
    b1 = ax_b.bar(x - w / 2, l2, w, color=colors, alpha=0.9, label="L2* disc.")
    ax2 = ax_b.twinx()
    b2 = ax2.bar(x + w / 2, nn, w, color=colors, alpha=0.55, hatch="///",
                 edgecolor="black", linewidth=0.4, label="NN_CV")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(labels, rotation=25, ha="right", fontsize=7)
    ax_b.set_ylabel("L2* discrepancy")
    ax2.set_ylabel("NN spacing CV")
    ax_b.set_title("(b) Coverage metrics", fontsize=10)
    ax_b.legend([b1, b2], ["L2* disc.", "NN_CV"], fontsize=7, loc="upper right")

    fig.tight_layout()
    return fig


def fig_si_coverage_vs_sample_size(
    method_curves: Dict[str, Tuple[np.ndarray, np.ndarray]],
    figsize: Tuple[float, float] = (5.5, 3.2),
) -> plt.Figure:
    """SI figure — coverage quality vs. ensemble size.

    Args:
        method_curves: {method_label: (n_array, metric_array)}. Typically
            produced by subsampling each method at sizes {50, 100, 200, 500}.
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize)
    palette = _default_colors_for_labels(list(method_curves.keys()))
    for (label, (ns, vals)), color in zip(method_curves.items(), palette):
        ax.plot(ns, vals, "o-", color=color, label=label, linewidth=1.4, markersize=4)
    ax.set_xlabel("Ensemble size")
    ax.set_ylabel("L2* discrepancy")
    ax.set_xscale("log")
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.set_title("Coverage quality vs. ensemble size")
    fig.tight_layout()
    return fig
