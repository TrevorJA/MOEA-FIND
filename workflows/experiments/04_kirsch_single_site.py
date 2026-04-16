"""Script 04 — Single-site Kirsch MOEA-FIND (manuscript Fig 5-6).

Couples Borg MOEA (or EpsNSGAII dev fallback) to the validated SynHydro
Kirsch-Nowak generator via KirschBorgWrapper. Supports both DV injection
modes (index, residual), SSI-based drought objectives, and bootstrap-
calibrated plausibility constraints.

Algorithm selection order:
    1. ``--algorithm`` CLI arg
    2. ``MOEA_FIND_ALGORITHM`` env var  (set by slurm scripts)
    3. Default: ``eps_nsga2`` (safe for local testing)

Constraint loading:
    If ``--constraints-json`` points to the output of
    ``scripts/diag_constraint_calibration.py``, calibrated tolerances are
    loaded. Otherwise constraints are disabled and a warning is printed.

Run locally (serial EpsNSGAII, unconstrained, 20 000 NFE):
    python scripts/04_kirsch_single_site.py --nfe 20000 --mode residual --plot

Run locally with placeholder constraints:
    python scripts/04_kirsch_single_site.py --nfe 20000 --mode residual \\
        --constraints-json outputs/diag_constraint_calibration/calibrated_tolerances.json

Run on HPC (MM Borg via MPI, constraints, 500 000 NFE):
    sbatch scripts/04_kirsch_single_site.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    run_experiment,
    make_variant_slug,
)
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.constraints import ConstraintConfig  # noqa: E402

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


def _load_constraints(
    json_path: Path | None,
    site_label: str,
    T_years: int,
) -> ConstraintConfig | None:
    """Load calibrated ConstraintConfig if JSON exists, else return None."""
    if json_path is None:
        return None
    json_path = Path(json_path)
    if not json_path.exists():
        print(f"[04] WARNING: constraints JSON not found at {json_path}; "
              f"running unconstrained.")
        return None
    try:
        cfg = ConstraintConfig.from_calibration_json(
            json_path, site_label=site_label, T_years=T_years,
        )
        print(f"[04] Loaded constraints from {json_path}")
        print(f"     annual_mean_tol={cfg.annual_mean_tol:.3f}, "
              f"annual_cv_tol={cfg.annual_cv_tol:.3f}, "
              f"lag1_ac_tol={cfg.lag1_ac_tol:.3f}, "
              f"non_drought_mean_tol={cfg.non_drought_mean_tol:.3f}, "
              f"seasonal_cycle_tol={cfg.seasonal_cycle_tol:.3f}")
        return cfg
    except Exception as exc:
        print(f"[04] WARNING: failed to load constraints: {exc}; "
              f"running unconstrained.")
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=20_000)
    p.add_argument("--n-years", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--mode", choices=["index", "residual"], default="residual")
    p.add_argument("--algorithm", default="eps_nsga2",
                   choices=["borg_mm", "borg_serial", "eps_nsga2"],
                   help="MOEA backend. Overridden by MOEA_FIND_ALGORITHM env var.")
    p.add_argument("--constraints-json", type=Path, default=None,
                   help="Path to calibrated_tolerances.json from Block A.")
    p.add_argument("--site-label", default="cannonsville")
    p.add_argument("--n-islands", type=int, default=1,
                   help="Number of islands for MM Borg.")
    p.add_argument("--checkpoint-freq", type=int, default=10000)
    p.add_argument("--old-checkpoint", type=Path, default=None,
                   help="Path to a checkpoint file to resume from.")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    # --- Constraints (load early — needed for variant slug) ---
    constraint_cfg = _load_constraints(
        args.constraints_json, args.site_label, args.n_years,
    )
    constrained = constraint_cfg is not None

    # --- Output directory keyed by variant slug ---
    slug = make_variant_slug(
        mode=args.mode,
        n_years=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        constrained=constrained,
    )
    out = args.output_dir / slug
    out.mkdir(parents=True, exist_ok=True)
    print(f"[04] variant: {slug}")

    # --- Data ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    # --- SSI characterisation ---
    objective_keys = ("mean_duration", "mean_avg_severity")
    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(
        monthly_1d, args.ssi
    )
    anti_ideal = compute_ssi_anti_ideal(hist_chars, objective_keys)
    print(f"[04] historical SSI-{args.ssi}: n={hist_chars['n_events']}, "
          f"dur={hist_chars['mean_duration']:.1f}mo, "
          f"severity={hist_chars['mean_avg_severity']:.3f}")
    print(f"[04] anti-ideal: {anti_ideal}")

    # --- Generator ---
    kirsch_gen = _build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(
        kirsch_gen, mode=args.mode, n_years_out=args.n_years,
    )

    # --- Config dump ---
    (out / "config.json").write_text(json.dumps({
        "script": "04_kirsch_single_site.py",
        "variant": slug,
        "algorithm": args.algorithm,
        "nfe": args.nfe,
        "n_years": args.n_years,
        "seed": args.seed,
        "ssi": args.ssi,
        "mode": args.mode,
        "constrained": constrained,
        "constraints_json": str(args.constraints_json) if args.constraints_json else None,
        "n_islands": args.n_islands,
    }, indent=2))

    # --- Run ---
    algo_kwargs = {}
    if args.algorithm in ("borg_mm", "borg_serial"):
        algo_kwargs["checkpoint_freq"] = args.checkpoint_freq
    if args.algorithm == "borg_mm":
        algo_kwargs["n_islands"] = args.n_islands
    if args.old_checkpoint:
        algo_kwargs["old_checkpoint"] = str(args.old_checkpoint)

    result = run_experiment(
        monthly_2d=monthly_2d,
        monthly_1d=monthly_1d,
        generator=generator,
        generator_name=f"Kirsch ({args.mode})",
        objective_keys=objective_keys,
        anti_ideal=anti_ideal,
        timescale=args.ssi,
        n_years_out=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        algorithm=args.algorithm,
        constraint_cfg=constraint_cfg,
        output_dir=out,
        **algo_kwargs,
    )

    (out / "results.json").write_text(
        json.dumps(result, indent=2, default=str)
    )
    print(f"[04] Pareto: {result.get('n_pareto', 0)} solutions")
    print(f"     wrote {out / 'results.json'}")

    # Save pareto.npz for downstream scripts (e.g., 10_plot_manuscript_figures)
    if result.get("n_pareto", 0) > 0:
        np.savez(
            out / "pareto.npz",
            dvs=np.array(result["pareto_dvs"]),
            objs=np.array(result["drought_metrics"]),
        )
        print(f"     wrote {out / 'pareto.npz'}")
    else:
        print(f"[04] WARNING: 0 Pareto solutions with {args.nfe} NFE "
              f"and {generator.n_dvs} DVs. Skipping pareto.npz and plots.")

    # --- Plots (saved to variant directory, not global figures/) ---
    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from src.plotting.drought_space import plot_scatter_with_marginals
            from src.plotting.trace_diagnostics import (
                plot_autocorrelation_comparison,
                plot_flow_duration_curve,
                plot_seasonal_cycle_comparison,
            )

            fig_dir = out / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)

            pareto_chars = result.get("pareto_chars", [])
            if pareto_chars:
                dm = np.array(result["drought_metrics"])
                hist_point = (
                    float(hist_chars["mean_duration"]),
                    float(hist_chars["mean_avg_severity"]),
                )
                fig_a = plot_scatter_with_marginals(
                    dm,
                    title=f"Kirsch ({args.mode}) Pareto",
                    historical_point=hist_point,
                    anti_ideal=anti_ideal,
                    objective_labels=(
                        "Mean duration (months)",
                        "Mean avg. severity",
                    ),
                )
                fig_a.savefig(fig_dir / "fig05_drought_space.pdf", dpi=300)
                plt.close(fig_a)
                print(f"[04] wrote {fig_dir / 'fig05_drought_space.pdf'}")

            pareto_traces_1d = result.get("pareto_traces_1d", [])
            pareto_traces_2d = result.get("pareto_traces_2d", [])
            if pareto_traces_1d:
                traces_1d = [np.array(t) for t in pareto_traces_1d]
                traces_2d = [np.array(t) for t in pareto_traces_2d]

                fig_acf, _ = plot_autocorrelation_comparison(
                    traces_1d, monthly_1d,
                )
                fig_acf.savefig(fig_dir / "fig06a_acf.pdf", dpi=300)
                plt.close(fig_acf)

                fig_fdc, _ = plot_flow_duration_curve(
                    traces_1d, monthly_1d,
                )
                fig_fdc.savefig(fig_dir / "fig06b_fdc.pdf", dpi=300)
                plt.close(fig_fdc)

                fig_sc, _ = plot_seasonal_cycle_comparison(
                    traces_2d, monthly_2d,
                )
                fig_sc.savefig(fig_dir / "fig06c_seasonal.pdf", dpi=300)
                plt.close(fig_sc)

                print(f"[04] wrote trace diagnostics to {fig_dir}")

        except ImportError as exc:
            print(f"[04] skipping plots (import error: {exc})")


if __name__ == "__main__":
    main()
