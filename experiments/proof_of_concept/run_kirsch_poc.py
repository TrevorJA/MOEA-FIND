"""MOEA-FIND Experiment: Kirsch nonparametric bootstrap with Borg wrapper.

Tests the KirschBorgWrapper class with both "index" and "residual" decision
variable modes using SSI-based drought objectives. The Kirsch method preserves
temporal autocorrelation, cross-year correlation (Dec-Jan), and seasonal
structure via Cholesky decomposition and normal-score transforms.

This experiment validates that the wrapper correctly maps Borg DVs to
indices/residuals and generates physically plausible synthetic traces.

Usage:
    python experiments/proof_of_concept/run_kirsch_poc.py --nfe 5000
    python experiments/proof_of_concept/run_kirsch_poc.py --mode residual --nfe 5000
"""

import argparse
import json
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.kirsch_wrapper import KirschBorgWrapper
from src.experiment_utils import (
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    run_experiment,
    plot_comparison,
)


def main():
    parser = argparse.ArgumentParser(
        description="Kirsch Borg wrapper proof of concept with SSI objectives"
    )
    parser.add_argument("--nfe", type=int, default=5000,
                        help="Number of function evaluations")
    parser.add_argument("--n-years", type=int, default=30,
                        help="Length of synthetic traces")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--ssi", type=int, default=3,
                        help="SSI accumulation period (1, 3, 6)")
    parser.add_argument("--mode", type=str, default="index",
                        choices=["index", "residual"],
                        help="KirschBorgWrapper mode")
    parser.add_argument("--compare", action="store_true",
                        help="Compare Kirsch index vs residual modes")
    args = parser.parse_args()

    output_dir = project_root / "outputs" / "proof_of_concept" / "kirsch"
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = project_root / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"

    # Load data
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    # Set up objectives and anti-ideal
    objective_keys = ("mean_duration", "mean_avg_severity")
    ssi_acc = args.ssi

    print(f"\n=== Kirsch Borg Wrapper PoC (SSI-{ssi_acc}) ===")

    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(
        monthly_1d, ssi_acc,
    )
    print(f"  Historical SSI-{ssi_acc}: {hist_chars['n_events']} events, "
          f"freq={hist_chars['frequency']:.1f}/dec, "
          f"dur={hist_chars['mean_duration']:.1f}mo, "
          f"severity={hist_chars['mean_avg_severity']:.3f}")

    anti_ideal = compute_ssi_anti_ideal(hist_chars, objective_keys)
    print(f"  Anti-ideal: dur={anti_ideal[0]:.1f}mo, severity={anti_ideal[1]:.3f}")

    if args.compare:
        # Compare Kirsch (index + residual) modes
        print(f"\n=== Comparison Mode ===")
        results_all = []

        # Kirsch generators (requires SynHydro KirschGenerator to be fitted first)
        try:
            import pandas as pd
            from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator

            # Create and fit Kirsch generator
            print(f"\n--- Kirsch Setup ---")
            kirsch_gen = KirschGenerator(generate_using_log_flow=True)

            # Convert monthly_2d to a DataFrame for Kirsch preprocessing
            dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
            daily_equiv = pd.DataFrame(
                {"flow_cfs": monthly_2d.flatten()},
                index=dates,
            )

            print(f"  Fitting Kirsch generator...")
            kirsch_gen.fit(daily_equiv)
            print(f"  Kirsch fitted: {kirsch_gen.n_historic_years} years, "
                  f"{kirsch_gen.n_sites} site(s)")

            for kirsch_mode in ["index", "residual"]:
                print(f"\n--- Kirsch ({kirsch_mode}) ---")
                generator = KirschBorgWrapper(
                    kirsch_gen, mode=kirsch_mode, n_years_out=args.n_years,
                )
                r = run_experiment(
                    monthly_2d, monthly_1d, generator,
                    f"Kirsch ({kirsch_mode})",
                    objective_keys, anti_ideal, ssi_acc,
                    args.n_years, args.nfe, args.seed,
                )
                results_all.append(r)

                fname = output_dir / f"results_ssi{ssi_acc}_kirsch_{kirsch_mode}.json"
                with open(fname, "w") as f:
                    json.dump(r, f, indent=2, default=str)

        except ImportError:
            print("\n  SynHydro not available; skipping Kirsch generators")

        # Generate comparison plots
        if len(results_all) >= 2:
            plot_comparison(
                results_all, hist_chars, anti_ideal, objective_keys, fig_dir,
            )

    else:
        # Single run of specified mode
        if args.mode in ("index", "residual"):
            try:
                import pandas as pd
                from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator

                print(f"\n--- Kirsch Setup ---")
                kirsch_gen = KirschGenerator(generate_using_log_flow=True)

                dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
                daily_equiv = pd.DataFrame(
                    {"flow_cfs": monthly_2d.flatten()},
                    index=dates,
                )

                print(f"  Fitting Kirsch generator...")
                kirsch_gen.fit(daily_equiv)

                print(f"\n--- Kirsch ({args.mode}) ---")
                generator = KirschBorgWrapper(
                    kirsch_gen, mode=args.mode, n_years_out=args.n_years,
                )

                r = run_experiment(
                    monthly_2d, monthly_1d, generator,
                    f"Kirsch ({args.mode})",
                    objective_keys, anti_ideal, ssi_acc,
                    args.n_years, args.nfe, args.seed,
                )

                fname = output_dir / f"results_ssi{ssi_acc}_kirsch_{args.mode}.json"
                with open(fname, "w") as f:
                    json.dump(r, f, indent=2, default=str)
                print(f"\nSaved: {fname}")

            except ImportError:
                print("\nERROR: SynHydro not available. Install via: pip install synhydro")
        else:
            print(f"ERROR: Unknown mode {args.mode!r}")


if __name__ == "__main__":
    main()
