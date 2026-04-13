"""Library-and-subsample baseline experiment for MOEA-FIND comparison.

Generates a large library of synthetic traces via SynHydro Kirsch,
characterizes each by SSI-based drought metrics, then subsamples using
LHS and Sobol designs in drought characteristic space.

This implements the conventional approach against which MOEA-FIND is compared.
See DD-09 in notes/design_decisions.md for context.

Usage:
    # Generate library (resource-intensive, run on HPC or patience)
    python experiments/kirsch_ensemble/run_library_baseline.py \\
        --n-traces 10000 --n-years 15 --ssi 3 --seed 42

    # Quick local test (small library)
    python experiments/kirsch_ensemble/run_library_baseline.py \\
        --n-traces 200 --n-years 15 --ssi 3 --seed 42

    # Subsample from existing library (fast, no generation)
    python experiments/kirsch_ensemble/run_library_baseline.py \\
        --load outputs/kirsch_ensemble/library/lib_10k \\
        --subsample-only --n-select 100
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.data import load_usgs_daily, daily_to_monthly
from src.library import LibraryGenerator
from src.objectives import (
    compute_ssi,
    compute_ssi_drought_characteristics,
    flows_to_series,
    make_ssi_calculator,
)
from src.analysis import coverage_metrics

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def prepare_data(cache_dir: Path):
    """Load and prepare historical USGS data.

    Returns:
        Tuple of (monthly_2d, monthly_1d) where:
            monthly_2d: shape (n_years, 12)
            monthly_1d: 1D array of all monthly flows
    """
    daily = load_usgs_daily(cache_dir=cache_dir)
    monthly = daily_to_monthly(daily)

    first_oct = monthly.index[monthly.index.month == 10][0]
    last_sep = monthly.index[monthly.index.month == 9][-1]
    monthly = monthly[first_oct:last_sep]

    n_years = len(monthly) // 12
    monthly_values = monthly.values[:n_years * 12]
    monthly_2d = monthly_values.reshape(n_years, 12)

    logger.info(
        "Historical: %d water years, mean=%.1f cfs",
        n_years, monthly_2d.mean(),
    )
    return monthly_2d, monthly_values


def fit_kirsch(monthly_2d: np.ndarray):
    """Create and fit a SynHydro KirschGenerator.

    Args:
        monthly_2d: Historical flows (n_years, 12).

    Returns:
        Fitted KirschGenerator instance.
    """
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator

    kirsch_gen = KirschGenerator(generate_using_log_flow=True)
    dates = pd.date_range(
        start="1950-10-01", periods=monthly_2d.size, freq="MS"
    )
    df = pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates)

    kirsch_gen.fit(df)
    logger.info(
        "Kirsch fitted: %d years, %d site(s)",
        kirsch_gen.n_historic_years, kirsch_gen.n_sites,
    )
    return kirsch_gen


def compute_historical_chars(monthly_1d, ssi_timescale, objective_keys):
    """Compute SSI drought characteristics for historical record.

    Returns:
        Dict of historical drought characteristics.
    """
    ssi_vals, calc = compute_ssi(monthly_1d, timescale=ssi_timescale)
    chars = compute_ssi_drought_characteristics(ssi_vals)
    logger.info(
        "Historical SSI-%d: %d events, freq=%.1f/dec, dur=%.1f mo",
        ssi_timescale, chars["n_events"],
        chars["frequency"], chars["mean_duration"],
    )
    return chars


def run_subsampling(
    lib: LibraryGenerator,
    objective_keys: tuple,
    n_select_values: list,
    seed: int = 42,
) -> list:
    """Run LHS and Sobol subsampling at multiple selection sizes.

    Args:
        lib: Characterized LibraryGenerator.
        objective_keys: Drought metrics for subsampling space.
        n_select_values: List of subsample sizes to try.
        seed: Random seed.

    Returns:
        List of result dicts (one per method x size combination).
    """
    results = []
    for n_sel in n_select_values:
        for method in ["lhs", "sobol"]:
            logger.info("Subsampling: %s, n=%d", method.upper(), n_sel)
            sub = lib.subsample(
                method=method,
                n_select=n_sel,
                objective_keys=objective_keys,
                seed=seed,
            )

            result = {
                "method": method,
                "n_requested": n_sel,
                "n_selected": sub["n_selected"],
                "coverage": sub["coverage"],
                "bounds": sub["bounds"],
                "objective_keys": list(objective_keys),
                "selected_points": sub["selected_points"].tolist(),
            }
            results.append(result)

            cov = sub["coverage"]
            l2 = cov.get("L2_star_discrepancy", 0)
            nncv = cov.get("nn_cv", 0)
            logger.info(
                "  %s n=%d: %d unique, L2*=%.4f, NN-CV=%.4f",
                method.upper(), n_sel, sub["n_selected"], l2, nncv,
            )

    return results


def plot_library_overview(lib: LibraryGenerator, objective_keys, fig_dir: Path):
    """Plot library drought characteristic distribution.

    Shows scatter of all library members and histograms of each metric.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chars = lib.characteristics
    k0, k1 = objective_keys

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Scatter of full library
    ax = axes[0]
    ax.scatter(chars[k0], chars[k1], s=3, alpha=0.3, c="#888888")
    ax.set_xlabel(k0.replace("_", " ").title())
    ax.set_ylabel(k1.replace("_", " ").title())
    ax.set_title(f"Full Library (n={len(chars)})")

    # Marginal histograms
    axes[1].hist(chars[k0], bins=50, color="#1f77b4", alpha=0.7, edgecolor="white")
    axes[1].set_xlabel(k0.replace("_", " ").title())
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"{k0} Distribution")

    axes[2].hist(chars[k1], bins=50, color="#ff7f0e", alpha=0.7, edgecolor="white")
    axes[2].set_xlabel(k1.replace("_", " ").title())
    axes[2].set_ylabel("Count")
    axes[2].set_title(f"{k1} Distribution")

    fig.suptitle(
        f"Kirsch Library: SSI-{lib.ssi_timescale} Drought Characteristics",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(fig_dir / "library_overview.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved library_overview.png")


def plot_subsample_comparison(
    lib: LibraryGenerator,
    subsample_results: list,
    objective_keys: tuple,
    fig_dir: Path,
):
    """Plot subsampled points overlaid on library.

    Shows how LHS and Sobol subsamples cover the feasible drought space.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chars = lib.characteristics
    k0, k1 = objective_keys

    # Group by n_requested
    sizes = sorted(set(r["n_requested"] for r in subsample_results))
    n_sizes = len(sizes)

    fig, axes = plt.subplots(n_sizes, 2, figsize=(12, 5 * n_sizes), squeeze=False)

    for row, n_sel in enumerate(sizes):
        for col, method in enumerate(["lhs", "sobol"]):
            ax = axes[row, col]

            # Background library
            ax.scatter(
                chars[k0], chars[k1],
                s=2, alpha=0.15, c="#cccccc", label="Library",
            )

            # Find matching result
            match = [
                r for r in subsample_results
                if r["method"] == method and r["n_requested"] == n_sel
            ]
            if match:
                pts = np.array(match[0]["selected_points"])
                cov = match[0]["coverage"]
                l2 = cov.get("L2_star_discrepancy", 0)
                n_actual = match[0]["n_selected"]

                ax.scatter(
                    pts[:, 0], pts[:, 1],
                    s=25, alpha=0.8,
                    c="#d62728" if method == "lhs" else "#1f77b4",
                    edgecolors="black", linewidths=0.3,
                    label=f"{method.upper()} (n={n_actual}, L2*={l2:.3f})",
                )

            ax.set_xlabel(k0.replace("_", " ").title())
            ax.set_ylabel(k1.replace("_", " ").title())
            ax.set_title(f"{method.upper()} subsample (target n={n_sel})")
            ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        f"Library Subsampling: SSI-{lib.ssi_timescale}",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(
        fig_dir / "library_subsample_comparison.png",
        dpi=200, bbox_inches="tight",
    )
    plt.close(fig)
    logger.info("Saved library_subsample_comparison.png")


def main():
    parser = argparse.ArgumentParser(
        description="Library-and-subsample baseline for MOEA-FIND comparison"
    )
    parser.add_argument(
        "--n-traces", type=int, default=200,
        help="Number of library traces to generate (default 200 for local testing)",
    )
    parser.add_argument("--n-years", type=int, default=15)
    parser.add_argument("--ssi", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--n-select", type=int, nargs="+", default=[50, 100, 200],
        help="Subsample sizes to evaluate",
    )
    parser.add_argument(
        "--load", type=str, default=None,
        help="Load existing library from this base path (skip generation)",
    )
    parser.add_argument(
        "--subsample-only", action="store_true",
        help="Only run subsampling (requires --load)",
    )
    args = parser.parse_args()

    output_dir = project_root / "outputs" / "kirsch_ensemble" / "library"
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = project_root / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    objective_keys = ("mean_duration", "mean_avg_severity")

    # Load historical data
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    hist_chars = compute_historical_chars(monthly_1d, args.ssi, objective_keys)

    # Fit Kirsch generator
    kirsch_gen = fit_kirsch(monthly_2d)

    # Create library generator
    lib = LibraryGenerator(
        kirsch_gen,
        n_years_out=args.n_years,
        ssi_timescale=args.ssi,
        reference_flows=monthly_1d,
    )

    lib_path = str(output_dir / f"lib_{args.n_traces // 1000}k")

    if args.load:
        # Load existing library
        lib.load_library(args.load)
        if lib.characteristics is None:
            logger.info("Characterizing loaded library...")
            lib.characterize(objective_keys=objective_keys)
    elif not args.subsample_only:
        # Generate and characterize
        t0 = time.time()
        lib.generate_library(n_traces=args.n_traces, seed=args.seed)
        gen_time = time.time() - t0
        logger.info("Generation time: %.1f s (%.1f ms/trace)", gen_time,
                     gen_time / args.n_traces * 1000)

        t0 = time.time()
        lib.characterize(objective_keys=objective_keys)
        char_time = time.time() - t0
        logger.info("Characterization time: %.1f s (%.1f ms/trace)", char_time,
                     char_time / args.n_traces * 1000)

        # Save library
        lib.save_library(lib_path)
    else:
        raise ValueError("--subsample-only requires --load")

    # Filter n_select values that are smaller than library size
    valid_n_select = [n for n in args.n_select if n <= len(lib.traces)]
    if not valid_n_select:
        logger.warning("All n_select values exceed library size %d", len(lib.traces))
        return

    # Run subsampling
    sub_results = run_subsampling(lib, objective_keys, valid_n_select, args.seed)

    # Save results
    results_path = output_dir / f"subsample_results_ssi{args.ssi}.json"
    with open(results_path, "w") as f:
        json.dump(sub_results, f, indent=2, default=str)
    logger.info("Saved results: %s", results_path)

    # Generate plots
    plot_library_overview(lib, objective_keys, fig_dir)
    plot_subsample_comparison(lib, sub_results, objective_keys, fig_dir)

    # Print summary
    print(f"\n=== Library Baseline Summary (SSI-{args.ssi}) ===")
    print(f"  Library size: {len(lib.traces)} traces")
    print(f"  Feasible region:")
    for key in objective_keys:
        vals = lib.characteristics[key]
        print(f"    {key}: [{vals.min():.2f}, {vals.max():.2f}]")

    print(f"\n  {'Method':<10} {'N target':>10} {'N unique':>10} "
          f"{'L2*':>10} {'NN-CV':>10}")
    print("  " + "-" * 55)
    for r in sub_results:
        cov = r["coverage"]
        print(
            f"  {r['method'].upper():<10} {r['n_requested']:>10} "
            f"{r['n_selected']:>10} "
            f"{cov.get('L2_star_discrepancy', 0):>10.4f} "
            f"{cov.get('nn_cv', 0):>10.4f}"
        )


if __name__ == "__main__":
    main()
