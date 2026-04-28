"""Script 13 — DV-uniformity ablation (manuscript Supporting Information).

Runs a single MOEA-FIND experiment under one of two mutually exclusive
constraint regimes so the two arms can be compared apples-to-apples:

    --arm hydrologic   : current 5-constraint plausibility formulation
                         (loaded from diag_constraint_calibration output).
    --arm dv_uniform   : single DV-space uniformity constraint with either
                         L2-star discrepancy or KS statistic
                         (loaded from diag_dv_uniformity_calibration output).

Both arms share the SAME objectives (DD-11: D_j + Manhattan to anti-ideal),
anti-ideal placement, DV length, seed, NFE, and algorithm — the only
difference is which constraint is active during the MOEA loop. Output is
written under ``outputs/exp13_dv_uniformity_ablation/<arm>/<slug>/`` so the
companion comparison script (``14_dv_uniformity_compare.py``) can collect
matched seed pairs.

Expected outcome: the DV-uniform arm will fail to find severe droughts,
confirming the falsification hypothesis documented in the plan file. That
result is the SI figure set.

Run locally (small NFE smoke test):
    python workflows/04_moea_find_single_site/dv_uniformity_ablation.py \\
        --arm dv_uniform --nfe 200 --seed 42 --algorithm eps_nsga2

Run on HPC (MM Borg, full NFE, array over arms x seeds):
    sbatch workflows/04_moea_find_single_site/slurm/dv_uniformity_ablation.slurm
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_config import DEFAULT_EXPERIMENT  # noqa: E402
from src.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.constraint_loaders import load_hydrologic_constraints, load_dv_uniformity_constraints  # noqa: E402
from src.experiment_utils import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    extract_pareto_maxes,
    run_experiment,
    make_variant_slug,
)
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.constraints import ConstraintConfig  # noqa: E402
from src.constraints_dv import DVUniformityConfig, VALID_STATISTICS  # noqa: E402

OUTPUT_SLUG = "exp13_dv_uniformity_ablation"
VALID_ARMS = ("hydrologic", "dv_uniform")


def main():
    cfg = DEFAULT_EXPERIMENT
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--arm", required=True, choices=VALID_ARMS,
                   help="Which constraint regime to run.")
    p.add_argument("--statistic", default=cfg.dv_uniformity_statistic,
                   choices=VALID_STATISTICS,
                   help="DV-uniformity statistic (ignored for hydrologic arm).")
    p.add_argument("--nfe", type=int, default=cfg.nfe)
    p.add_argument("--n-years", type=int, default=cfg.n_years_out)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--ssi", type=int, default=cfg.ssi_timescale,
                   choices=[1, 3, 6, 12])
    p.add_argument("--mode", choices=["index", "residual"], default=cfg.dv_mode)
    p.add_argument("--algorithm", default=cfg.algorithm,
                   choices=["borg_mm", "borg_serial", "eps_nsga2"])
    p.add_argument("--constraints-json", type=Path, default=cfg.constraints_json,
                   help="Hydrologic calibrated_tolerances.json (hydrologic arm).")
    p.add_argument("--dv-uniformity-json", type=Path,
                   default=cfg.dv_uniformity_json,
                   help="DV-uniformity calibrated_dv_tolerances.json "
                        "(dv_uniform arm).")
    p.add_argument("--site-label", default=cfg.site_label)
    p.add_argument("--anti-ideal-reference", type=Path,
                   default=cfg.anti_ideal_reference_json,
                   help="Prior results.json whose Pareto max drives D* "
                        "for non-cyclic objectives. Cyclic metrics remain "
                        "at 12*headroom.")
    p.add_argument("--n-islands", type=int, default=cfg.n_islands)
    p.add_argument("--checkpoint-freq", type=int, default=cfg.checkpoint_freq)
    p.add_argument("--plot", action="store_true",
                   help="Render per-arm diagnostic figures locally "
                        "(comparison figures are built by script 14).")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    args = p.parse_args()

    # --- Load the arm-appropriate constraint ---
    hydrologic_cfg: ConstraintConfig | None = None
    dv_cfg: DVUniformityConfig | None = None
    if args.arm == "hydrologic":
        hydrologic_cfg = load_hydrologic_constraints(
            args.constraints_json, args.site_label, args.n_years,
        )
    else:  # dv_uniform
        dv_cfg = load_dv_uniformity_constraints(
            args.dv_uniformity_json, args.site_label, args.n_years,
            args.statistic,
        )
    constrained = hydrologic_cfg is not None or dv_cfg is not None

    # --- Output directory: arm first, then the usual variant slug ---
    extra = {"cm": args.arm}
    if args.arm == "dv_uniform":
        extra["st"] = args.statistic
    slug = make_variant_slug(
        mode=args.mode,
        n_years=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        constrained=constrained,
        extra=extra,
    )
    out = args.output_dir / args.arm / slug
    out.mkdir(parents=True, exist_ok=True)
    print(f"[13] arm={args.arm} slug={slug}")

    # --- Data ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    # --- SSI characterisation ---
    objective_keys = cfg.objective_keys
    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(
        monthly_1d, args.ssi
    )
    feasible_maxes = None
    if args.anti_ideal_reference is not None:
        ref = Path(args.anti_ideal_reference)
        if ref.exists():
            feasible_maxes = extract_pareto_maxes(ref, objective_keys)
            print(f"[13] anti-ideal reference from {ref}")
            print(f"     Pareto maxes: {feasible_maxes}")
        else:
            print(f"[13] WARNING: --anti-ideal-reference {ref} not found; "
                  f"falling back to historical-max D*.")

    anti_ideal = compute_ssi_anti_ideal(
        hist_chars, objective_keys, headroom=cfg.anti_ideal_headroom,
        feasible_maxes=feasible_maxes,
    )
    print(f"[13] anti-ideal: {anti_ideal}")

    # --- Generator ---
    kirsch_gen = build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(
        kirsch_gen, mode=args.mode, n_years_out=args.n_years,
    )

    # --- Sanity: n_dvs vs calibration ---
    if dv_cfg is not None and generator.n_dvs != dv_cfg.n_dvs:
        print(f"[13] WARNING: generator n_dvs={generator.n_dvs} differs from "
              f"calibrated dv_cfg.n_dvs={dv_cfg.n_dvs}. "
              f"Rebuild the calibration with matching --T / --dv-mode, "
              f"or compute_dv_constraint will raise at evaluation time.")

    # --- Config dump ---
    (out / "config.json").write_text(json.dumps({
        "script": "13_dv_uniformity_ablation.py",
        "variant": slug,
        "arm": args.arm,
        "statistic": args.statistic if args.arm == "dv_uniform" else None,
        "algorithm": args.algorithm,
        "nfe": args.nfe,
        "n_years": args.n_years,
        "seed": args.seed,
        "ssi": args.ssi,
        "mode": args.mode,
        "constrained": constrained,
        "constraints_json": (str(args.constraints_json)
                             if args.constraints_json else None),
        "dv_uniformity_json": (str(args.dv_uniformity_json)
                               if args.dv_uniformity_json else None),
        "n_islands": args.n_islands,
    }, indent=2))

    # --- Run ---
    algo_kwargs = {}
    if args.algorithm in ("borg_mm", "borg_serial"):
        algo_kwargs["checkpoint_freq"] = args.checkpoint_freq
    if args.algorithm == "borg_mm":
        algo_kwargs["n_islands"] = args.n_islands

    result = run_experiment(
        monthly_2d=monthly_2d,
        monthly_1d=monthly_1d,
        generator=generator,
        generator_name=f"Kirsch ({args.mode}) / {args.arm}",
        objective_keys=objective_keys,
        anti_ideal=anti_ideal,
        timescale=args.ssi,
        n_years_out=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        algorithm=args.algorithm,
        constraint_cfg=hydrologic_cfg,
        dv_constraint_cfg=dv_cfg,
        output_dir=out,
        **algo_kwargs,
    )

    n_pareto = result.get("n_pareto", 0)
    print(f"[13] arm={args.arm} Pareto: {n_pareto} solutions")

    # Only the MM Borg master rank has Pareto solutions; workers must not
    # overwrite master output files.
    if n_pareto > 0:
        (out / "results.json").write_text(
            json.dumps(result, indent=2, default=str)
        )
        print(f"     wrote {out / 'results.json'}")
        np.savez(
            out / "pareto.npz",
            dvs=np.array(result["pareto_dvs"]),
            objs=np.array(result["drought_metrics"]),
        )
        print(f"     wrote {out / 'pareto.npz'}")
    else:
        print(f"[13] arm={args.arm}: 0 Pareto solutions; skipping output "
              f"writes (either MM Borg worker rank or empty archive).")

    # Per-arm diagnostic plots are optional here; the comparison figures
    # for the SI are produced by workflows/04_moea_find_single_site/dv_uniformity_compare.py.
    if args.plot and n_pareto > 0:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from src.plotting.drought_space import plot_scatter_with_marginals
            from src.historical_blocks import compute_historical_block_chars

            fig_dir = out / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)

            dm = np.array(result["drought_metrics"])
            hist_point = (
                float(hist_chars["mean_duration"]),
                float(hist_chars["mean_avg_severity"]),
            )
            hist_block_chars = compute_historical_block_chars(
                monthly_1d, T_years=args.n_years, ssi_calc=ssi_calc,
                objective_keys=objective_keys, stride=1,
            )
            fig_a = plot_scatter_with_marginals(
                dm[:, :2],
                title=f"{args.arm} Pareto vs historical blocks",
                historical_point=hist_point,
                anti_ideal=anti_ideal[:2],
                historical_cloud=hist_block_chars[:, :2],
                objective_labels=(
                    "Mean duration (months)",
                    "Mean avg. severity",
                ),
            )
            fig_a.savefig(fig_dir / "drought_space.pdf", dpi=200)
            plt.close(fig_a)
        except ImportError as exc:
            print(f"[13] skipping per-arm plots (import error: {exc})")


if __name__ == "__main__":
    main()
