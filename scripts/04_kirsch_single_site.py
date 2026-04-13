"""Script 04 — Single-site Kirsch MOEA-FIND (manuscript §6.1-6.2, Fig 5-6).

Couples MM Borg MOEA to the validated SynHydro Kirsch-Nowak generator via
KirschBorgWrapper. Supports both DV injection modes (index, residual) and
SSI-based drought objectives (mean_duration, mean_avg_severity) plus the
Manhattan-norm auxiliary objective.

Run locally (serial Borg/platypus fallback):
    python scripts/04_kirsch_single_site.py --nfe 20000 --mode residual --plot

Run on HPC (MM Borg via MPI):
    sbatch scripts/04_kirsch_single_site.slurm

Outputs under outputs/exp04_kirsch_single_site/<mode>/:
    - results.json  (Pareto objectives, drought characteristics, timing)
    - pareto.npz    (DVs, objectives)
    - config.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    run_experiment,
)
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402

OUTPUT_SLUG = "exp04_kirsch_single_site"


def _build_kirsch_generator(monthly_2d: np.ndarray):
    """Fit SynHydro's KirschGenerator on historical monthly flows."""
    import pandas as pd
    from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator

    gen = KirschGenerator(generate_using_log_flow=True)
    dates = pd.date_range(start="1950-10-01", periods=monthly_2d.size, freq="MS")
    daily_equiv = pd.DataFrame({"flow_cfs": monthly_2d.flatten()}, index=dates)
    gen.fit(daily_equiv)
    return gen


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=20_000)
    p.add_argument("--n-years", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--mode", choices=["index", "residual"], default="residual")
    p.add_argument("--constrained", action="store_true",
                   help="Add lag-1 ACF and non-drought annual stats constraints (Exp 2.2).")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    out = args.output_dir / args.mode
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps({
        "script": "04_kirsch_single_site.py",
        "manuscript_section": (
            "§6.2 Constrained Pareto (Fig 6)" if args.constrained
            else "§6.1 Unconstrained Pareto (Fig 5)"
        ),
        "nfe": args.nfe, "n_years": args.n_years, "seed": args.seed,
        "ssi": args.ssi, "mode": args.mode, "constrained": args.constrained,
    }, indent=2))

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    objective_keys = ("mean_duration", "mean_avg_severity")
    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(monthly_1d, args.ssi)
    anti_ideal = compute_ssi_anti_ideal(hist_chars, objective_keys)
    print(f"[04] historical SSI-{args.ssi}: n={hist_chars['n_events']}, "
          f"dur={hist_chars['mean_duration']:.1f}mo, "
          f"severity={hist_chars['mean_avg_severity']:.3f}")
    print(f"[04] anti-ideal: {anti_ideal}")

    kirsch_gen = _build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(kirsch_gen, mode=args.mode, n_years_out=args.n_years)

    result = run_experiment(
        monthly_2d, monthly_1d, generator,
        f"Kirsch ({args.mode})",
        objective_keys, anti_ideal, args.ssi,
        args.n_years, args.nfe, args.seed,
    )

    (out / "results.json").write_text(json.dumps(result, indent=2, default=str))
    print(f"[04] Pareto: {len(result.get('pareto_objs', []))} solutions")
    print(f"     wrote {out / 'results.json'}")

    if args.plot:
        from src.plotting.drought_space import plot_scatter_with_marginals
        from src.plotting.trace_diagnostics import (
            plot_autocorrelation_comparison, plot_flow_duration_curve,
            plot_seasonal_cycle_comparison,
        )
        import matplotlib.pyplot as plt

        pareto_chars = np.asarray(result["pareto_drought_chars"])
        hist_point = (float(hist_chars["mean_duration"]),
                      float(hist_chars["mean_avg_severity"]))
        fig_a = plot_scatter_with_marginals(
            pareto_chars, title=f"§6 Kirsch ({args.mode}) Pareto",
            historical_point=hist_point, anti_ideal=anti_ideal,
            objective_labels=("Mean duration (months)", "Mean avg. severity"),
        )
        fig_a.savefig(PROJECT_ROOT / "figures" / f"fig05_kirsch_{args.mode}.pdf",
                      dpi=300)
        plt.close(fig_a)

        synth_traces = [np.asarray(t) for t in result.get("synthetic_traces", [])]
        if synth_traces:
            fig_b, _ = plot_autocorrelation_comparison(
                synth_traces, monthly_1d.flatten(), max_lag=24,
            )
            fig_b.savefig(PROJECT_ROOT / "figures" / f"fig06a_acf_{args.mode}.pdf",
                          dpi=300)
            plt.close(fig_b)
            fig_c, _ = plot_flow_duration_curve(synth_traces, monthly_1d.flatten())
            fig_c.savefig(PROJECT_ROOT / "figures" / f"fig06b_fdc_{args.mode}.pdf",
                          dpi=300)
            plt.close(fig_c)
            synth_2d = [t.reshape(-1, 12) for t in synth_traces
                        if t.size % 12 == 0]
            if synth_2d:
                fig_d, _ = plot_seasonal_cycle_comparison(synth_2d, monthly_2d)
                fig_d.savefig(PROJECT_ROOT / "figures"
                              / f"fig06c_seasonal_{args.mode}.pdf", dpi=300)
                plt.close(fig_d)


if __name__ == "__main__":
    main()
