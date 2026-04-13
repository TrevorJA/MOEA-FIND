"""Diagnostic figures for Kirsch ensemble experiments.

Generates publication-quality figures from existing result JSON files.
Does NOT regenerate traces or rerun optimization.

Figures produced:
  1. Pareto front comparison (index vs residual) in drought space
  2. Drought characteristic correlation matrix (all metrics)
  3. Coverage quality comparison (L2*, NN-CV, with LHS/Sobol baselines)
  4. Hyperplane verification (objective sum distributions)
  5. Pareto front with marginal histograms
  6. Multi-metric pairwise scatter (all drought chars)
  7. Pareto density and spacing analysis

Usage:
    python experiments/kirsch_ensemble/plot_diagnostics.py
    python experiments/kirsch_ensemble/plot_diagnostics.py --results-dir outputs/kirsch_ensemble
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.analysis import coverage_metrics, generate_lhs_samples, generate_sobol_samples


def load_results(results_dir: Path) -> dict:
    """Load all Kirsch result JSON files from a directory.

    Args:
        results_dir: Path to directory containing result JSONs.

    Returns:
        Dict mapping mode name to result dict.
    """
    results = {}
    for f in sorted(results_dir.glob("results_ssi*_kirsch_*.json")):
        data = json.load(open(f))
        results[data["mode"]] = data
    return results


def extract_pareto_chars_matrix(result: dict) -> tuple:
    """Extract a matrix of all drought characteristics from pareto_chars.

    Args:
        result: Single experiment result dict.

    Returns:
        Tuple of (matrix of shape (n_pareto, n_metrics), list of metric names).
    """
    chars_list = result["pareto_chars"]
    if not chars_list:
        return np.array([]), []

    # Use consistent metric ordering
    metric_keys = [
        "frequency", "mean_duration", "mean_magnitude",
        "mean_severity", "mean_avg_severity",
        "max_duration", "max_magnitude", "worst_severity",
    ]
    # Filter to metrics that exist
    available = [k for k in metric_keys if k in chars_list[0]]

    matrix = np.array([[c.get(k, 0.0) for k in available] for c in chars_list])
    return matrix, available


def fig1_pareto_comparison(results: dict, fig_dir: Path):
    """Pareto front scatter: index vs residual in objective space.

    Shows both modes overlaid with historical point and anti-ideal.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    colors = {"Kirsch (index)": "#d62728", "Kirsch (residual)": "#9467bd"}
    markers = {"Kirsch (index)": "o", "Kirsch (residual)": "s"}

    for ax_idx, (mode, r) in enumerate(results.items()):
        ax = axes[ax_idx]
        metrics = np.array(r["drought_metrics"])
        obj_keys = r["objective_keys"]
        anti_ideal = np.array(r["anti_ideal"])
        n = r["n_pareto"]
        cov = r["coverage"]
        l2 = cov.get("L2_star_discrepancy", cov.get("L2_star", 0))

        ax.scatter(
            metrics[:, 0], metrics[:, 1],
            s=20, alpha=0.6, c=colors.get(mode, "gray"),
            marker=markers.get(mode, "o"), edgecolors="none",
        )
        ax.scatter(
            anti_ideal[0], anti_ideal[1],
            marker="X", s=150, c="black", zorder=5, label="Anti-ideal",
        )
        ax.set_xlabel("Mean Duration (months)", fontsize=11)
        ax.set_ylabel("Mean Avg Severity (SSI units)", fontsize=11)
        ax.set_title(f"{mode}\n(n={n}, L2*={l2:.3f})", fontsize=11)
        ax.legend(fontsize=9)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    fig.suptitle(
        f"SSI-{list(results.values())[0]['ssi_timescale']} "
        "Drought Space: Kirsch Index vs Residual",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_pareto_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig1_pareto_comparison.png")


def fig2_correlation_matrix(results: dict, fig_dir: Path):
    """Correlation matrix of all drought characteristics across Pareto solutions.

    Shows how drought metrics co-vary. Highlights the duration-magnitude
    correlation that motivated using severity instead.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_modes = len(results)
    fig, axes = plt.subplots(1, n_modes, figsize=(7 * n_modes, 6))
    if n_modes == 1:
        axes = [axes]

    for ax, (mode, r) in zip(axes, results.items()):
        matrix, keys = extract_pareto_chars_matrix(r)
        if matrix.size == 0:
            continue

        # Compute correlation
        # Filter out zero-variance columns
        stds = matrix.std(axis=0)
        valid = stds > 1e-10
        matrix_valid = matrix[:, valid]
        keys_valid = [k for k, v in zip(keys, valid) if v]

        corr = np.corrcoef(matrix_valid.T)

        im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
        ax.set_xticks(range(len(keys_valid)))
        ax.set_yticks(range(len(keys_valid)))

        # Shorter labels
        short = {
            "frequency": "freq", "mean_duration": "dur",
            "mean_magnitude": "mag", "mean_severity": "sev",
            "mean_avg_severity": "avg_sev", "max_duration": "max_dur",
            "max_magnitude": "max_mag", "worst_severity": "worst_sev",
        }
        labels = [short.get(k, k) for k in keys_valid]
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        ax.set_yticklabels(labels, fontsize=9)

        # Annotate cells
        for i in range(len(keys_valid)):
            for j in range(len(keys_valid)):
                color = "white" if abs(corr[i, j]) > 0.6 else "black"
                ax.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center",
                        fontsize=7, color=color)

        ax.set_title(f"{mode} (n={r['n_pareto']})", fontsize=11)

    fig.colorbar(im, ax=axes, shrink=0.8, label="Pearson r")
    fig.suptitle("Drought Characteristic Correlations (Pareto Solutions)", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_correlation_matrix.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig2_correlation_matrix.png")


def fig3_coverage_comparison(results: dict, fig_dir: Path):
    """Coverage quality: MOEA-FIND vs LHS and Sobol baselines.

    For each mode, computes L2* and NN-CV, then generates matched-size
    LHS and Sobol samples for comparison.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modes = []
    l2_moea = []
    l2_lhs = []
    l2_sobol = []
    nncv_moea = []
    nncv_lhs = []
    nncv_sobol = []

    for mode, r in results.items():
        metrics = np.array(r["drought_metrics"])
        n = len(metrics)
        d = metrics.shape[1]
        anti_ideal = np.array(r["anti_ideal"])
        lb = np.zeros(d)
        ub = anti_ideal.copy()

        cov_moea = coverage_metrics(metrics, lb, ub)

        # Generate baselines with same n
        lhs_pts = generate_lhs_samples(n, d, lb, ub, seed=42)
        sobol_pts = generate_sobol_samples(n, d, lb, ub, seed=42)[:n]
        cov_lhs = coverage_metrics(lhs_pts, lb, ub)
        cov_sobol = coverage_metrics(sobol_pts, lb, ub)

        modes.append(mode)
        l2_moea.append(cov_moea["L2_star_discrepancy"])
        l2_lhs.append(cov_lhs["L2_star_discrepancy"])
        l2_sobol.append(cov_sobol["L2_star_discrepancy"])
        nncv_moea.append(cov_moea["nn_cv"])
        nncv_lhs.append(cov_lhs["nn_cv"])
        nncv_sobol.append(cov_sobol["nn_cv"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    x = np.arange(len(modes))
    w = 0.25

    # L2* discrepancy (lower is better)
    ax1.bar(x - w, l2_moea, w, label="MOEA-FIND", color="#d62728")
    ax1.bar(x, l2_lhs, w, label="LHS baseline", color="#2ca02c")
    ax1.bar(x + w, l2_sobol, w, label="Sobol baseline", color="#1f77b4")
    ax1.set_xticks(x)
    ax1.set_xticklabels([m.replace("Kirsch ", "") for m in modes])
    ax1.set_ylabel("L2* Discrepancy (lower = more uniform)")
    ax1.set_title("Space-Filling Quality")
    ax1.legend(fontsize=9)

    # NN-CV (lower is better = more regular spacing)
    ax2.bar(x - w, nncv_moea, w, label="MOEA-FIND", color="#d62728")
    ax2.bar(x, nncv_lhs, w, label="LHS baseline", color="#2ca02c")
    ax2.bar(x + w, nncv_sobol, w, label="Sobol baseline", color="#1f77b4")
    ax2.set_xticks(x)
    ax2.set_xticklabels([m.replace("Kirsch ", "") for m in modes])
    ax2.set_ylabel("NN-CV (lower = more regular spacing)")
    ax2.set_title("Spacing Regularity")
    ax2.legend(fontsize=9)

    fig.suptitle("Coverage Quality: MOEA-FIND vs Space-Filling Baselines", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_coverage_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig3_coverage_comparison.png")


def fig4_hyperplane_verification(results: dict, fig_dir: Path):
    """Verify Manhattan norm forces solutions onto hyperplane.

    Histogram of sum(objectives) for each Pareto solution. Should be
    a delta function at sum(anti_ideal).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_modes = len(results)
    fig, axes = plt.subplots(1, n_modes, figsize=(6 * n_modes, 4))
    if n_modes == 1:
        axes = [axes]

    for ax, (mode, r) in zip(axes, results.items()):
        metrics = np.array(r["drought_metrics"])
        anti_ideal = np.array(r["anti_ideal"])
        obj_keys = r["objective_keys"]

        # Reconstruct Manhattan norm for each solution
        manhattan = np.sum(np.abs(metrics - anti_ideal), axis=1)
        obj_sums = np.sum(metrics, axis=1) + manhattan
        expected = np.sum(anti_ideal)

        ax.hist(obj_sums, bins=50, color=("#d62728" if "index" in mode else "#9467bd"),
                alpha=0.8, edgecolor="white")
        ax.axvline(expected, color="black", linestyle="--", linewidth=2,
                   label=f"Expected = {expected:.2f}")
        ax.set_xlabel("Sum of all objectives", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"{mode}\nstd = {r['hyperplane']['actual_std']:.6f}", fontsize=11)
        ax.legend(fontsize=9)

    fig.suptitle("Hyperplane Verification: J1 + J2 + J_Manhattan = constant", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_hyperplane.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig4_hyperplane.png")


def fig5_marginal_scatter(results: dict, fig_dir: Path):
    """Joint scatter with marginal histograms for each mode.

    Publication-style figure showing the Pareto front distribution
    with marginal density on each axis.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    for mode, r in results.items():
        metrics = np.array(r["drought_metrics"])
        obj_keys = r["objective_keys"]
        n = r["n_pareto"]

        fig = plt.figure(figsize=(7, 7))
        gs = GridSpec(4, 4, hspace=0.05, wspace=0.05)

        ax_main = fig.add_subplot(gs[1:4, 0:3])
        ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
        ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

        color = "#d62728" if "index" in mode else "#9467bd"

        ax_main.scatter(metrics[:, 0], metrics[:, 1], s=15, alpha=0.6,
                        c=color, edgecolors="none")
        ax_main.set_xlabel("Mean Duration (months)", fontsize=11)
        ax_main.set_ylabel("Mean Avg Severity (SSI units)", fontsize=11)

        ax_top.hist(metrics[:, 0], bins=30, color=color, alpha=0.7, edgecolor="white")
        ax_top.set_ylabel("Count", fontsize=9)
        plt.setp(ax_top.get_xticklabels(), visible=False)

        ax_right.hist(metrics[:, 1], bins=30, color=color, alpha=0.7,
                      edgecolor="white", orientation="horizontal")
        ax_right.set_xlabel("Count", fontsize=9)
        plt.setp(ax_right.get_yticklabels(), visible=False)

        safe_mode = mode.replace(" ", "_").replace("(", "").replace(")", "")
        fig.suptitle(f"{mode} (n={n})", fontsize=13)
        fig.savefig(fig_dir / f"fig5_marginal_{safe_mode}.png",
                    dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved fig5_marginal_{safe_mode}.png")


def fig6_pairwise_chars(results: dict, fig_dir: Path):
    """Pairwise scatter of all drought characteristics.

    Shows how the Pareto ensemble spans multiple drought dimensions,
    not just the two optimized objectives.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Pick one mode with more solutions
    best_mode = max(results.keys(), key=lambda m: results[m]["n_pareto"])
    r = results[best_mode]
    matrix, keys = extract_pareto_chars_matrix(r)

    if matrix.size == 0:
        return

    # Select interesting subset
    show_keys = ["frequency", "mean_duration", "mean_avg_severity",
                 "max_duration", "worst_severity"]
    show_idx = [i for i, k in enumerate(keys) if k in show_keys]
    show_labels = [keys[i] for i in show_idx]
    matrix_sub = matrix[:, show_idx]
    d = len(show_idx)

    short = {
        "frequency": "Frequency\n(events/dec)",
        "mean_duration": "Mean Duration\n(months)",
        "mean_avg_severity": "Mean Avg\nSeverity",
        "max_duration": "Max Duration\n(months)",
        "worst_severity": "Worst\nSeverity",
    }

    fig, axes = plt.subplots(d, d, figsize=(2.5 * d, 2.5 * d))
    color = "#d62728" if "index" in best_mode else "#9467bd"

    for i in range(d):
        for j in range(d):
            ax = axes[i, j]
            if i == j:
                ax.hist(matrix_sub[:, i], bins=20, color=color,
                        alpha=0.7, edgecolor="white")
            elif i > j:
                ax.scatter(matrix_sub[:, j], matrix_sub[:, i],
                           s=8, alpha=0.4, c=color, edgecolors="none")
            else:
                # Upper triangle: correlation value
                valid_mask = ~(np.isnan(matrix_sub[:, i]) | np.isnan(matrix_sub[:, j]))
                if valid_mask.sum() > 2:
                    corr = np.corrcoef(matrix_sub[valid_mask, i],
                                       matrix_sub[valid_mask, j])[0, 1]
                    ax.text(0.5, 0.5, f"r={corr:.2f}", ha="center", va="center",
                            fontsize=12, transform=ax.transAxes,
                            fontweight="bold" if abs(corr) > 0.7 else "normal",
                            color="red" if abs(corr) > 0.7 else "black")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)

            if i == d - 1:
                ax.set_xlabel(short.get(show_labels[j], show_labels[j]), fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(short.get(show_labels[i], show_labels[i]), fontsize=8)
            else:
                ax.set_yticklabels([])

    fig.suptitle(f"Pairwise Drought Characteristics: {best_mode}", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig6_pairwise_chars.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig6_pairwise_chars.png")


def fig7_spacing_analysis(results: dict, fig_dir: Path):
    """Nearest-neighbor distance distribution for Pareto solutions.

    Shows how uniformly the Pareto front tiles the drought space.
    Includes kernel density estimate and comparison across modes.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.spatial import KDTree

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    colors = {"Kirsch (index)": "#d62728", "Kirsch (residual)": "#9467bd"}

    for mode, r in results.items():
        metrics = np.array(r["drought_metrics"])
        anti_ideal = np.array(r["anti_ideal"])

        # Normalize to [0, 1]
        lb = metrics.min(axis=0)
        ub = metrics.max(axis=0)
        rng = ub - lb
        rng[rng == 0] = 1
        normed = (metrics - lb) / rng

        tree = KDTree(normed)
        nn_dists, _ = tree.query(normed, k=2)
        nn_dists = nn_dists[:, 1]  # Nearest neighbor (not self)

        c = colors.get(mode, "gray")
        ax1.hist(nn_dists, bins=30, alpha=0.6, color=c, label=mode, edgecolor="white")

        # Cumulative distribution
        sorted_d = np.sort(nn_dists)
        cdf = np.arange(1, len(sorted_d) + 1) / len(sorted_d)
        ax2.plot(sorted_d, cdf, color=c, linewidth=2, label=mode)

    ax1.set_xlabel("Nearest-Neighbor Distance (normalized)", fontsize=11)
    ax1.set_ylabel("Count", fontsize=11)
    ax1.set_title("NN Distance Distribution", fontsize=11)
    ax1.legend(fontsize=9)

    ax2.set_xlabel("Nearest-Neighbor Distance (normalized)", fontsize=11)
    ax2.set_ylabel("Cumulative Probability", fontsize=11)
    ax2.set_title("NN Distance CDF", fontsize=11)
    ax2.legend(fontsize=9)

    fig.suptitle("Pareto Front Spacing Analysis", fontsize=13)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig7_spacing_analysis.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig7_spacing_analysis.png")


def print_summary_table(results: dict):
    """Print a formatted summary table of all results."""
    print("\n  === Kirsch Ensemble Diagnostics Summary ===")
    print(f"  {'Mode':<22} {'N':>5} {'Duration':>14} {'Severity':>14} "
          f"{'L2*':>8} {'NN-CV':>8} {'Hyperplane':>12}")
    print("  " + "-" * 90)
    for mode, r in results.items():
        rng = r["ranges"]
        cov = r["coverage"]
        obj_keys = r["objective_keys"]
        k0, k1 = obj_keys[0], obj_keys[1]
        l2 = cov.get("L2_star_discrepancy", 0)
        nncv = cov.get("nn_cv", 0)
        hp_std = r["hyperplane"]["actual_std"]
        print(
            f"  {mode:<22} {r['n_pareto']:>5} "
            f"{rng[k0]['min']:.1f}-{rng[k0]['max']:.1f} mo  "
            f"{rng[k1]['min']:.2f}-{rng[k1]['max']:.2f}     "
            f"{l2:>8.4f} {nncv:>8.4f} "
            f"{'PASS' if hp_std < 0.01 else f'std={hp_std:.3f}':>12}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Generate diagnostic figures from Kirsch ensemble results"
    )
    parser.add_argument(
        "--results-dir", type=str,
        default=str(project_root / "outputs" / "kirsch_ensemble"),
        help="Directory containing result JSON files",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading results from: {results_dir}")
    results = load_results(results_dir)

    if not results:
        print("ERROR: No result files found. Run run_kirsch_experiment.py first.")
        return

    print(f"Found {len(results)} result sets: {list(results.keys())}")

    print_summary_table(results)

    print("\nGenerating figures...")
    fig1_pareto_comparison(results, fig_dir)
    fig2_correlation_matrix(results, fig_dir)
    fig3_coverage_comparison(results, fig_dir)
    fig4_hyperplane_verification(results, fig_dir)
    fig5_marginal_scatter(results, fig_dir)
    fig6_pairwise_chars(results, fig_dir)
    fig7_spacing_analysis(results, fig_dir)

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == "__main__":
    main()
