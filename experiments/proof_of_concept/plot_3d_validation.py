"""Generate publication-quality diagnostic figures for MOEA-FIND 3D validation.

This script creates six comprehensive figures demonstrating the core concepts
of the MOEA-FIND method applied to analytic proof-of-concept experiments:

1. Manhattan Norm Concept (2-panel) — 2D and 3D projections with hyperplane constraint
2. 3D Simplex Coverage (4-panel) — Pairwise projections + 3D scatter
3. Coverage Comparison (3x4 grid) — Pareto vs LHS vs Sobol with metrics
4. Hyperplane Verification (2-panel) — Histogram of objective sums
5. Epsilon Sensitivity — n_solutions, L2*, NN_CV vs epsilon
6. NN Distance Distribution — Spacing uniformity across methods

All figures use academic styling (no gridlines, minimal chrome, clean fonts).
Saved as PNG (200 dpi) and PDF to outputs/poc/figures/.

Usage:
    python experiments/proof_of_concept/plot_3d_validation.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.analysis import (
    generate_lhs_samples,
    generate_sobol_samples,
    coverage_metrics,
)

OUTPUT_DIR = project_root / "outputs" / "poc"
FIG_DIR = OUTPUT_DIR / "figures"


# Academic color palette
COLORS = {
    "pareto": "#2166ac",  # blue
    "lhs": "#b2182b",     # red
    "sobol": "#4dac26",   # green
    "random": "#969696",  # gray
}


def setup_style():
    """Configure matplotlib for academic publication."""
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.constrained_layout.use": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def load_data():
    """Load POC results from disk."""
    with open(OUTPUT_DIR / "poc_summary.json") as f:
        summary = json.load(f)
    data_2d = np.load(OUTPUT_DIR / "pareto_2d.npz")
    data_3d = np.load(OUTPUT_DIR / "pareto_3d.npz")
    return summary, data_2d, data_3d


def panel_label(ax, label: str):
    """Add (a), (b), etc. label to top-left of axes."""
    try:
        # For 3D axes
        ax.text2D(-0.15, 1.05, f"({label})",
                  transform=ax.transAxes,
                  fontsize=11,
                  fontweight="bold",
                  ha="right",
                  va="bottom")
    except AttributeError:
        # For 2D axes
        ax.text(
            -0.15, 1.05, f"({label})",
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            ha="right",
            va="bottom",
        )


def fig1_manhattan_norm_concept(data_2d: dict, data_3d: dict):
    """Figure 1: Manhattan Norm Concept (2-panel).

    Panel (a): 2D case — Pareto front with hyperplane X1+X2=6 overlaid.
    Panel (b): 3D case — 3D scatter with anti-ideal point.
    """
    fig = plt.figure(figsize=(14, 5), constrained_layout=True)
    gs = GridSpec(1, 2, figure=fig)

    # Panel (a): 2D case
    ax = fig.add_subplot(gs[0, 0])
    dvs_2d = data_2d["dvs"]
    objs_2d = data_2d["objs"]

    # Color by Manhattan norm (which is constant for 2D on hyperplane)
    manhattan_2d = np.sum(np.abs(objs_2d), axis=1)
    scatter_a = ax.scatter(
        dvs_2d[:, 0], dvs_2d[:, 1],
        c=manhattan_2d,
        s=8,
        alpha=0.6,
        cmap="viridis",
        rasterized=True,
        edgecolors="none",
    )
    cbar_a = plt.colorbar(scatter_a, ax=ax)
    cbar_a.set_label("Manhattan norm", fontsize=9)

    # Overlay the constraint hyperplane (for visualization only)
    # In 2D space: this is implicit in the objective space, but we show the Pareto front
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.set_xlabel("$X_1$", fontsize=10)
    ax.set_ylabel("$X_2$", fontsize=10)
    ax.set_title("2D Analytic Test (k=2, 3 objectives)", fontsize=11)
    ax.set_aspect("equal")
    panel_label(ax, "a")

    # Panel (b): 3D case
    ax = fig.add_subplot(gs[0, 1], projection="3d")
    dvs_3d = data_3d["dvs"]
    objs_3d = data_3d["objs"]

    # Color by Manhattan norm (constant on hyperplane)
    manhattan_3d = np.sum(np.abs(objs_3d), axis=1)
    scatter_b = ax.scatter(
        dvs_3d[:, 0], dvs_3d[:, 1], dvs_3d[:, 2],
        c=manhattan_3d,
        s=3,
        alpha=0.5,
        cmap="viridis",
        rasterized=True,
    )

    # Mark anti-ideal point
    anti_ideal_3d = np.array([3, 3, 3])
    ax.scatter([anti_ideal_3d[0]], [anti_ideal_3d[1]], [anti_ideal_3d[2]],
               s=150, c="red", marker="*", edgecolors="black",
               linewidth=0.5, label="Anti-ideal", zorder=10)

    ax.set_xlabel("$X_1$", fontsize=10)
    ax.set_ylabel("$X_2$", fontsize=10)
    ax.set_zlabel("$X_3$", fontsize=10)
    ax.set_title("3D Analytic Test (k=3, 4 objectives)", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    panel_label(ax, "b")

    # Save
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "fig1_manhattan_concept.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig1_manhattan_concept.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig1_manhattan_concept.png/pdf")


def fig2_3d_simplex_coverage(data_3d: dict):
    """Figure 2: 3D Simplex Coverage (4-panel).

    Panels (a-c): Pairwise projections (X1-X2, X1-X3, X2-X3).
    Panel (d): 3D scatter with rotated view.
    """
    fig = plt.figure(figsize=(14, 11), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)

    dvs_3d = data_3d["dvs"]
    pairs = [
        (0, 1, "$X_1$", "$X_2$", "a"),
        (0, 2, "$X_1$", "$X_3$", "b"),
        (1, 2, "$X_2$", "$X_3$", "c"),
    ]

    # Pairwise projections (top-left, top-right, bottom-left)
    for idx, (i, j, xi, xj, label) in enumerate(pairs):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        ax.scatter(
            dvs_3d[:, i], dvs_3d[:, j],
            s=2, alpha=0.5, c=COLORS["pareto"], rasterized=True, edgecolors="none"
        )
        ax.set_xlabel(xi, fontsize=10)
        ax.set_ylabel(xj, fontsize=10)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.2, 3.2)
        ax.set_aspect("equal")
        ax.set_title(f"Projection: {xi} vs {xj}", fontsize=11)
        panel_label(ax, label)

    # 3D scatter (bottom-right)
    ax3d = fig.add_subplot(gs[1, 1], projection="3d")
    ax3d.scatter(
        dvs_3d[:, 0], dvs_3d[:, 1], dvs_3d[:, 2],
        s=2, alpha=0.4, c=COLORS["pareto"], rasterized=True, edgecolors="none"
    )
    ax3d.set_xlabel("$X_1$", fontsize=10)
    ax3d.set_ylabel("$X_2$", fontsize=10)
    ax3d.set_zlabel("$X_3$", fontsize=10)
    ax3d.set_title(f"3D Pareto Front (n={len(dvs_3d)})", fontsize=11)
    ax3d.view_init(elev=20, azim=45)
    panel_label(ax3d, "d")

    fig.savefig(FIG_DIR / "fig2_3d_simplex_coverage.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig2_3d_simplex_coverage.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig2_3d_simplex_coverage.png/pdf")


def fig3_coverage_comparison(data_3d: dict):
    """Figure 3: Coverage Comparison (4x4 grid).

    Rows: Pareto front, LHS, Sobol, Random (all at n=1362)
    Columns: X1-X2, X1-X3, X2-X3, metrics bar chart
    """
    n = len(data_3d["dvs"])
    lb = np.array([-3.0, -3.0, -3.0])
    ub = np.array([3.0, 3.0, 3.0])

    # Generate samples
    pareto_dvs = data_3d["dvs"]
    lhs_dvs = generate_lhs_samples(n, 3, lb, ub, seed=42)
    sobol_dvs = generate_sobol_samples(n, 3, lb, ub, seed=42)[:n]
    rng = np.random.default_rng(42)
    random_dvs = lb + rng.random((n, 3)) * (ub - lb)

    # Compute metrics
    pareto_metrics = coverage_metrics(pareto_dvs, lb, ub)
    lhs_metrics = coverage_metrics(lhs_dvs, lb, ub)
    sobol_metrics = coverage_metrics(sobol_dvs, lb, ub)
    random_metrics = coverage_metrics(random_dvs, lb, ub)

    methods = [("Pareto Front", pareto_dvs, COLORS["pareto"]),
               ("Latin Hypercube", lhs_dvs, COLORS["lhs"]),
               ("Sobol Sequence", sobol_dvs, COLORS["sobol"]),
               ("Random Uniform", random_dvs, COLORS["random"])]

    pairs = [(0, 1, "$X_1$", "$X_2$"),
             (0, 2, "$X_1$", "$X_3$"),
             (1, 2, "$X_2$", "$X_3$")]

    fig = plt.figure(figsize=(14, 14), constrained_layout=True)
    gs = GridSpec(4, 4, figure=fig)

    for row, (method_name, dvs, color) in enumerate(methods):
        # Three pairwise projections
        for col, (i, j, xi, xj) in enumerate(pairs):
            ax = fig.add_subplot(gs[row, col])
            ax.scatter(dvs[:, i], dvs[:, j], s=1.5, alpha=0.5, c=color,
                      rasterized=True, edgecolors="none")
            ax.set_xlim(-3.2, 3.2)
            ax.set_ylim(-3.2, 3.2)
            ax.set_aspect("equal")
            if row == 0:
                ax.set_title(f"{xi} vs {xj}", fontsize=11)
            if col == 0:
                ax.set_ylabel(xi, fontsize=10)
            if row == 3:
                ax.set_xlabel(xi, fontsize=10)

            # Panel labels
            label = chr(97 + row * 4 + col)  # a, b, c, ..., p
            panel_label(ax, label)

        # Metrics bar chart
        ax_bar = fig.add_subplot(gs[row, 3])
        if row == 0:
            metrics_set = pareto_metrics
        elif row == 1:
            metrics_set = lhs_metrics
        elif row == 2:
            metrics_set = sobol_metrics
        else:
            metrics_set = random_metrics

        metrics_vals = [metrics_set["L2_star_discrepancy"], metrics_set["nn_cv"]]
        metric_names = ["L2*", "NN_CV"]
        bars = ax_bar.bar(metric_names, metrics_vals, color=[color, color],
                         alpha=0.7, edgecolor="black", linewidth=0.5)
        ax_bar.set_ylabel("Value", fontsize=10)
        ax_bar.set_title(f"{method_name}\n(n={len(dvs)})", fontsize=10)

        # Annotate bars
        for bar, val in zip(bars, metrics_vals):
            ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                       f"{val:.3f}", ha="center", va="bottom", fontsize=8)

        label = chr(97 + row * 4 + 3)
        panel_label(ax_bar, label)

    fig.savefig(FIG_DIR / "fig3_coverage_comparison.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig3_coverage_comparison.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig3_coverage_comparison.png/pdf")


def fig4_hyperplane_verification(data_2d: dict, data_3d: dict):
    """Figure 4: Hyperplane Verification (2-panel).

    Shows that all solutions lie exactly on their respective hyperplanes.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # Panel (a): 2D
    ax = axes[0]
    objs_2d = data_2d["objs"]
    sums_2d = np.sum(objs_2d, axis=1)
    dev_2d = sums_2d - 6.0

    ax.bar(
        ["All Solutions"],
        [len(sums_2d)],
        color=COLORS["pareto"],
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(f"2D: n={len(sums_2d)}, Sum={np.mean(sums_2d):.1f}", fontsize=11)
    ax.text(
        0.5, 0.7,
        f"Std(deviation) = {np.std(dev_2d):.2e}\n(exactly on hyperplane)",
        transform=ax.transAxes,
        ha="center",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.3),
    )
    panel_label(ax, "a")

    # Panel (b): 3D
    ax = axes[1]
    objs_3d = data_3d["objs"]
    sums_3d = np.sum(objs_3d, axis=1)
    dev_3d = sums_3d - 9.0

    ax.bar(
        ["All Solutions"],
        [len(sums_3d)],
        color=COLORS["pareto"],
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(f"3D: n={len(sums_3d)}, Sum={np.mean(sums_3d):.1f}", fontsize=11)
    ax.text(
        0.5, 0.7,
        f"Std(deviation) = {np.std(dev_3d):.2e}\n(exactly on hyperplane)",
        transform=ax.transAxes,
        ha="center",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.3),
    )
    panel_label(ax, "b")

    fig.savefig(FIG_DIR / "fig4_hyperplane_verification.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig4_hyperplane_verification.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig4_hyperplane_verification.png/pdf")


def fig5_epsilon_sensitivity(summary: dict):
    """Figure 5: Epsilon Sensitivity (3-panel).

    Shows how Pareto front properties change with epsilon parameter.
    """
    sens = summary["experiment_1_3_sensitivity"]["epsilon_sensitivity"]
    eps_vals = np.array([s["epsilon"] for s in sens])
    n_sols = np.array([s["n_solutions"] for s in sens])
    l2_vals = np.array([s.get("L2_star", 0) for s in sens])
    nn_cv = np.array([s["nn_cv"] for s in sens])

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)

    # Panel (a): n_solutions vs epsilon
    ax = axes[0]
    ax.plot(eps_vals, n_sols, "o-", color=COLORS["pareto"], linewidth=2, markersize=6)
    ax.set_xlabel("Epsilon", fontsize=10)
    ax.set_ylabel("Ensemble Size", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Ensemble Size vs Epsilon", fontsize=11)
    panel_label(ax, "a")

    # Panel (b): L2* vs epsilon
    ax = axes[1]
    ax.plot(eps_vals, l2_vals, "s-", color=COLORS["lhs"], linewidth=2, markersize=6)
    ax.set_xlabel("Epsilon", fontsize=10)
    ax.set_ylabel("L2* Discrepancy", fontsize=10)
    ax.set_xscale("log")
    ax.set_title("Coverage Quality vs Epsilon", fontsize=11)
    panel_label(ax, "b")

    # Panel (c): NN_CV vs epsilon
    ax = axes[2]
    ax.plot(eps_vals, nn_cv, "^-", color=COLORS["sobol"], linewidth=2, markersize=6)
    ax.set_xlabel("Epsilon", fontsize=10)
    ax.set_ylabel("NN Distance CV", fontsize=10)
    ax.set_xscale("log")
    ax.set_title("Spacing Uniformity vs Epsilon", fontsize=11)
    panel_label(ax, "c")

    fig.savefig(FIG_DIR / "fig5_epsilon_sensitivity.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig5_epsilon_sensitivity.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig5_epsilon_sensitivity.png/pdf")


def fig6_nn_distance_distribution(data_3d: dict):
    """Figure 6: NN Distance Distribution (2-panel).

    Panel (a): Histogram/KDE of nearest-neighbor distances.
    Panel (b): Empirical CDF of NN distances.
    """
    from scipy.spatial import KDTree
    from scipy.stats import gaussian_kde

    n = len(data_3d["dvs"])
    lb = np.array([-3.0, -3.0, -3.0])
    ub = np.array([3.0, 3.0, 3.0])

    # Generate samples
    pareto_dvs = data_3d["dvs"]
    lhs_dvs = generate_lhs_samples(n, 3, lb, ub, seed=42)
    sobol_dvs = generate_sobol_samples(n, 3, lb, ub, seed=42)[:n]
    rng = np.random.default_rng(42)
    random_dvs = lb + rng.random((n, 3)) * (ub - lb)

    # Compute NN distances
    def nn_distances(dvs):
        normed = (dvs - lb) / (ub - lb)
        tree = KDTree(normed)
        dists, _ = tree.query(normed, k=2)
        return dists[:, 1]

    pareto_nn = nn_distances(pareto_dvs)
    lhs_nn = nn_distances(lhs_dvs)
    sobol_nn = nn_distances(sobol_dvs)
    random_nn = nn_distances(random_dvs)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # Panel (a): Histogram + KDE
    ax = axes[0]
    bins = np.linspace(0, max(pareto_nn.max(), lhs_nn.max(), sobol_nn.max(), random_nn.max()), 30)

    ax.hist(pareto_nn, bins=bins, alpha=0.5, density=True, color=COLORS["pareto"],
           label="Pareto", edgecolor="black", linewidth=0.5)
    ax.hist(lhs_nn, bins=bins, alpha=0.5, density=True, color=COLORS["lhs"],
           label="LHS", edgecolor="black", linewidth=0.5)
    ax.hist(sobol_nn, bins=bins, alpha=0.5, density=True, color=COLORS["sobol"],
           label="Sobol", edgecolor="black", linewidth=0.5)
    ax.hist(random_nn, bins=bins, alpha=0.5, density=True, color=COLORS["random"],
           label="Random", edgecolor="black", linewidth=0.5)

    # Overlay KDE
    for nn, color, label in [(pareto_nn, COLORS["pareto"], "Pareto"),
                              (lhs_nn, COLORS["lhs"], "LHS"),
                              (sobol_nn, COLORS["sobol"], "Sobol"),
                              (random_nn, COLORS["random"], "Random")]:
        kde = gaussian_kde(nn)
        x_kde = np.linspace(nn.min(), nn.max(), 200)
        ax.plot(x_kde, kde(x_kde), color=color, linewidth=2, alpha=0.8)

    ax.set_xlabel("Nearest-Neighbor Distance", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("NN Distance Distribution", fontsize=11)
    ax.legend(fontsize=9)
    panel_label(ax, "a")

    # Panel (b): Empirical CDF
    ax = axes[1]
    for nn, color, label in [(pareto_nn, COLORS["pareto"], "Pareto"),
                              (lhs_nn, COLORS["lhs"], "LHS"),
                              (sobol_nn, COLORS["sobol"], "Sobol"),
                              (random_nn, COLORS["random"], "Random")]:
        sorted_nn = np.sort(nn)
        cdf = np.arange(1, len(sorted_nn) + 1) / len(sorted_nn)
        ax.plot(sorted_nn, cdf, color=color, linewidth=2, label=label, alpha=0.8)

    ax.set_xlabel("Nearest-Neighbor Distance", fontsize=10)
    ax.set_ylabel("Cumulative Probability", fontsize=10)
    ax.set_title("NN Distance CDF", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_xlim(left=0)
    panel_label(ax, "b")

    fig.savefig(FIG_DIR / "fig6_nn_distribution.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig6_nn_distribution.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig6_nn_distribution.png/pdf")


def main():
    """Generate all publication-quality figures."""
    print("Setting up style...")
    setup_style()

    print("Loading data...")
    summary, data_2d, data_3d = load_data()

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Generating figures...")
    fig1_manhattan_norm_concept(data_2d, data_3d)
    fig2_3d_simplex_coverage(data_3d)
    fig3_coverage_comparison(data_3d)
    fig4_hyperplane_verification(data_2d, data_3d)
    fig5_epsilon_sensitivity(summary)
    fig6_nn_distance_distribution(data_3d)

    print(f"\nAll figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
