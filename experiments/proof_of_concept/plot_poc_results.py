"""Generate figures for MOEA-FIND Proof-of-Concept experiments.

Figures:
1. 2D Pareto front projected to (X1, X2) vs LHS vs Sobol
2. 3D Pareto front projected to (X1, X2, X3) — 3 pairwise projections + 3D scatter
3. Epsilon sensitivity: n_solutions and discrepancy vs epsilon
4. Hyperplane verification: objective sums histogram

Usage:
    python experiments/proof_of_concept/plot_poc_results.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.analysis import generate_lhs_samples, generate_sobol_samples, coverage_metrics

OUTPUT_DIR = project_root / "outputs" / "poc"
FIG_DIR = project_root / "outputs" / "poc" / "figures"


def load_results():
    """Load POC results."""
    with open(OUTPUT_DIR / "poc_summary.json") as f:
        summary = json.load(f)
    data_2d = np.load(OUTPUT_DIR / "pareto_2d.npz")
    data_3d = np.load(OUTPUT_DIR / "pareto_3d.npz")
    return summary, data_2d, data_3d


def fig1_2d_coverage_comparison(pareto_dvs: np.ndarray, seed: int = 42):
    """Figure 1: 2D Pareto front vs LHS vs Sobol vs Random."""
    n = len(pareto_dvs)
    lb = np.array([-3.0, -3.0])
    ub = np.array([3.0, 3.0])

    lhs = generate_lhs_samples(n, 2, lb, ub, seed=seed)
    sobol = generate_sobol_samples(n, 2, lb, ub, seed=seed)[:n]
    rng = np.random.default_rng(seed)
    rand = lb + rng.random((n, 2)) * (ub - lb)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharex=True, sharey=True)

    datasets = [
        ("MOEA-FIND\n(Pareto Front)", pareto_dvs, "C0"),
        ("Latin Hypercube", lhs, "C1"),
        ("Sobol Sequence", sobol, "C2"),
        ("Random Uniform", rand, "C3"),
    ]

    for ax, (title, data, color) in zip(axes, datasets):
        ax.scatter(data[:, 0], data[:, 1], s=1.5, alpha=0.5, c=color, rasterized=True)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.2, 3.2)
        ax.set_aspect("equal")

        # Compute and display discrepancy
        normed = (data - lb) / (ub - lb)
        from scipy.stats.qmc import discrepancy
        d = discrepancy(normed, method="L2-star")
        ax.text(0.05, 0.95, f"L2*={d:.5f}\nn={len(data)}",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    axes[0].set_xlabel("$X_1$")
    axes[0].set_ylabel("$X_2$")
    for ax in axes[1:]:
        ax.set_xlabel("$X_1$")

    fig.suptitle("Experiment 1.1: Coverage Comparison (2D)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig1_2d_coverage.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved fig1_2d_coverage.png")


def fig2_3d_coverage(pareto_dvs: np.ndarray):
    """Figure 2: 3D Pareto front — pairwise projections and 3D scatter."""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    pairs = [(0, 1, "$X_1$", "$X_2$"), (0, 2, "$X_1$", "$X_3$"), (1, 2, "$X_2$", "$X_3$")]

    # Pairwise 2D projections (top row)
    for idx, (i, j, xi, xj) in enumerate(pairs):
        ax = fig.add_subplot(gs[0, idx])
        ax.scatter(pareto_dvs[:, i], pareto_dvs[:, j], s=1, alpha=0.4, c="C0", rasterized=True)
        ax.set_xlabel(xi)
        ax.set_ylabel(xj)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.2, 3.2)
        ax.set_aspect("equal")
        ax.set_title(f"Projection: {xi} vs {xj}", fontsize=10)

    # 3D scatter (bottom left, spanning 2 columns)
    ax3d = fig.add_subplot(gs[1, 0:2], projection="3d")
    ax3d.scatter(
        pareto_dvs[:, 0], pareto_dvs[:, 1], pareto_dvs[:, 2],
        s=1, alpha=0.3, c="C0", rasterized=True,
    )
    ax3d.set_xlabel("$X_1$")
    ax3d.set_ylabel("$X_2$")
    ax3d.set_zlabel("$X_3$")
    ax3d.set_title(f"3D Pareto Front (n={len(pareto_dvs)})", fontsize=10)

    # Coverage comparison bar chart (bottom right)
    n = len(pareto_dvs)
    lb = np.array([-3.0, -3.0, -3.0])
    ub = np.array([3.0, 3.0, 3.0])
    lhs = generate_lhs_samples(n, 3, lb, ub)
    sobol = generate_sobol_samples(n, 3, lb, ub)[:n]

    pareto_m = coverage_metrics(pareto_dvs, lb, ub)
    lhs_m = coverage_metrics(lhs, lb, ub)
    sobol_m = coverage_metrics(sobol, lb, ub)

    ax_bar = fig.add_subplot(gs[1, 2])
    methods = ["Pareto", "LHS", "Sobol"]
    l2_vals = [pareto_m["L2_star_discrepancy"], lhs_m["L2_star_discrepancy"], sobol_m["L2_star_discrepancy"]]
    colors = ["C0", "C1", "C2"]
    bars = ax_bar.bar(methods, l2_vals, color=colors, alpha=0.7)
    ax_bar.set_ylabel("L2* Discrepancy")
    ax_bar.set_title("3D Coverage Quality", fontsize=10)
    for bar, val in zip(bars, l2_vals):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{val:.5f}", ha="center", va="bottom", fontsize=8)

    fig.suptitle("Experiment 1.2: 3D Analytic Test (k=3, 4 objectives)", fontsize=12, y=1.01)
    fig.savefig(FIG_DIR / "fig2_3d_coverage.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved fig2_3d_coverage.png")


def fig3_epsilon_sensitivity(summary: dict):
    """Figure 3: Epsilon sensitivity analysis."""
    sens = summary["experiment_1_3_sensitivity"]["epsilon_sensitivity"]
    eps_vals = [s["epsilon"] for s in sens]
    n_sols = [s["n_solutions"] for s in sens]
    l2_vals = [s.get("L2_star_discrepancy", s.get("L2_star", 0)) for s in sens]
    nn_cv = [s["nn_cv"] for s in sens]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))

    ax1.plot(eps_vals, n_sols, "o-", color="C0")
    ax1.set_xlabel("Epsilon")
    ax1.set_ylabel("Pareto Solutions")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title("Ensemble Size vs Epsilon")
    ax1.grid(True, alpha=0.3)

    ax2.plot(eps_vals, l2_vals, "s-", color="C1")
    ax2.set_xlabel("Epsilon")
    ax2.set_ylabel("L2* Discrepancy")
    ax2.set_xscale("log")
    ax2.set_title("Coverage Quality vs Epsilon")
    ax2.grid(True, alpha=0.3)

    ax3.plot(eps_vals, nn_cv, "^-", color="C2")
    ax3.set_xlabel("Epsilon")
    ax3.set_ylabel("NN Distance CV")
    ax3.set_xscale("log")
    ax3.set_title("Spacing Uniformity vs Epsilon")
    ax3.grid(True, alpha=0.3)

    fig.suptitle("Experiment 1.3: Epsilon Sensitivity (2D)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig3_epsilon_sensitivity.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved fig3_epsilon_sensitivity.png")


def fig4_hyperplane_verification(data_2d, data_3d):
    """Figure 4: Verify that objective sums are constant (hyperplane property)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    objs_2d = data_2d["objs"]
    sums_2d = np.sum(objs_2d, axis=1)
    dev_2d = sums_2d - 6.0
    if np.std(dev_2d) < 1e-10:
        ax1.bar(["All solutions"], [len(sums_2d)], color="C0", alpha=0.7)
        ax1.text(0.5, 0.7, f"Sum = {np.mean(sums_2d):.6f}\n(exactly on hyperplane)",
                 transform=ax1.transAxes, ha="center", fontsize=9)
    else:
        ax1.hist(dev_2d, bins=50, color="C0", alpha=0.7, edgecolor="black", linewidth=0.5)
    ax1.set_xlabel("Deviation from expected sum")
    ax1.set_ylabel("Count")
    ax1.set_title(f"2D: n={len(sums_2d)}, expected $\\Sigma$=6.0")

    objs_3d = data_3d["objs"]
    sums_3d = np.sum(objs_3d, axis=1)
    dev_3d = sums_3d - 9.0
    if np.std(dev_3d) < 1e-10:
        ax2.bar(["All solutions"], [len(sums_3d)], color="C1", alpha=0.7)
        ax2.text(0.5, 0.7, f"Sum = {np.mean(sums_3d):.6f}\n(exactly on hyperplane)",
                 transform=ax2.transAxes, ha="center", fontsize=9)
    else:
        ax2.hist(dev_3d, bins=50, color="C1", alpha=0.7, edgecolor="black", linewidth=0.5)
    ax2.set_xlabel("Deviation from expected sum")
    ax2.set_ylabel("Count")
    ax2.set_title(f"3D: n={len(sums_3d)}, expected $\\Sigma$=9.0")

    fig.suptitle("Hyperplane Verification: $\\sum J_i = \\sum D^*_i$", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig4_hyperplane_check.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved fig4_hyperplane_check.png")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading results...")
    summary, data_2d, data_3d = load_results()

    pareto_2d = data_2d["dvs"]
    pareto_3d = data_3d["dvs"]

    print("Generating figures...")
    fig1_2d_coverage_comparison(pareto_2d)
    fig2_3d_coverage(pareto_3d)
    fig3_epsilon_sensitivity(summary)
    fig4_hyperplane_verification(data_2d, data_3d)
    print("Done.")


if __name__ == "__main__":
    main()
