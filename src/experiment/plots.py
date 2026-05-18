"""Cross-generator comparison figure for MOEA-FIND experiments.

Split out of the former ``src.experiment_utils`` god module.
"""

from pathlib import Path

import numpy as np


def plot_comparison(
    results_list: list,
    hist_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys,
    fig_dir: Path,
) -> None:
    """Generate comparison figures for multiple generators.

    Creates scatter plots of Pareto fronts in drought space, overlaying
    historical point and anti-ideal point.

    Args:
        results_list: List of results dicts from run_experiment.
        hist_chars: Historical drought characteristics dict.
        anti_ideal: Anti-ideal point (k-dimensional).
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.metrics.drought_metrics.DroughtMetric` instances.
        fig_dir: Output directory for figures.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.metrics.drought_metrics import resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    if len(metric_set) < 2:
        raise ValueError("plot_comparison requires at least two metrics")

    m0, m1 = metric_set[0], metric_set[1]

    fig_dir.mkdir(parents=True, exist_ok=True)

    color_map = {
        "Kirsch (index)": "#d62728",
        "Kirsch (residual)": "#9467bd",
    }

    n_modes = len(results_list)
    fig, axes = plt.subplots(1, n_modes, figsize=(6 * n_modes, 5), squeeze=False)
    axes = axes[0]

    for ax, r in zip(axes, results_list):
        metrics = np.array(r["drought_metrics"])
        mode = r["mode"]
        c = color_map.get(mode, "gray")

        ax.scatter(
            metrics[:, 0], metrics[:, 1],
            s=15, alpha=0.7, c=c,
            label=f"{mode} (n={r['n_pareto']})",
        )
        ax.scatter(
            m0.extract(hist_chars), m1.extract(hist_chars),
            marker="*", s=200, c="black", zorder=5, label="Historical",
        )
        ax.scatter(
            anti_ideal[0], anti_ideal[1],
            marker="x", s=200, c="red", zorder=5, label="Anti-ideal D*",
        )
        ax.set_xlabel(f"{m0.label} ({m0.units})")
        ax.set_ylabel(f"{m1.label} ({m1.units})")
        ax.set_title(f"{mode} (n={r['n_pareto']})")
        ax.legend(fontsize=8)

    ssi_acc = results_list[0].get("ssi_timescale", "?")
    fig.suptitle(f"Kirsch: SSI-{ssi_acc} Drought Space (index vs residual)", fontsize=13)
    fig.tight_layout()

    fname = f"kirsch_poc_ssi{ssi_acc}_coverage.png"
    fig.savefig(fig_dir / fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname}")

    # Summary table
    print(f"\n  === SSI-{ssi_acc} Kirsch PoC Summary ===")
    print(f"  {'Generator':<25} {'N':>5} "
          f"{m0.label[:16]:>16} {m1.label[:16]:>16} {'L2*':>8}")
    for r in results_list:
        rng = r["ranges"]
        cov = r["coverage"]
        l2 = cov.get("L2_star_discrepancy", cov.get("L2_star", 0))
        r0 = rng.get(m0.name, {"min": 0.0, "max": 0.0})
        r1 = rng.get(m1.name, {"min": 0.0, "max": 0.0})
        print(f"  {r['mode']:<25} {r['n_pareto']:>5} "
              f"{r0['min']:>7.2f}-{r0['max']:<7.2f} "
              f"{r1['min']:>7.2f}-{r1['max']:<7.2f} "
              f"{l2:>8.4f}")
