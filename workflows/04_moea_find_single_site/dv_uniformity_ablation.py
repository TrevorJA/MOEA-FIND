"""Stage 04 / dv_uniformity_ablation -- SI-G arm (compute only).

Runs a single MOEA-FIND experiment under one constraint regime so the
arms can be compared apples-to-apples by ``dv_uniformity_compare.py``.

Arms:
    --arm hydrologic   : 5-statistic plausibility formulation
    --arm dv_uniform   : DV-space uniformity constraint with --statistic
                         (l2_star or ad).

Both arms share objectives (DD-11), anti-ideal placement, DV length,
NFE, and algorithm. Output written under
``outputs/04_moea_find_single_site/dv_uniformity_ablation/<arm>/<slug>/``.
Figures are produced by the paired
``workflows/99_supporting_info_figures/dv_uniformity_compare.py``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment.config import DEFAULT_EXPERIMENT  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.optimization.constraint_loaders import load_hydrologic_constraints, load_dv_uniformity_constraints  # noqa: E402
from src.experiment import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    extract_pareto_maxes,
    run_experiment,
)
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.optimization.constraints import ConstraintConfig  # noqa: E402
from src.optimization.constraints_dv import DVUniformityConfig, VALID_STATISTICS  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import moea_slug  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "dv_uniformity_ablation"
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
                   choices=["borg_mm", "borg_serial"])
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
                        "for non-cyclic objectives.")
    p.add_argument("--n-islands", type=int, default=cfg.n_islands)
    p.add_argument("--checkpoint-freq", type=int, default=cfg.checkpoint_freq)
    args = p.parse_args()

    # --- Load the arm-appropriate constraint ---
    hydrologic_cfg: ConstraintConfig | None = None
    dv_cfg: DVUniformityConfig | None = None
    if args.arm == "hydrologic":
        hydrologic_cfg = load_hydrologic_constraints(
            args.constraints_json, args.site_label, args.n_years,
        )
    else:
        dv_cfg = load_dv_uniformity_constraints(
            args.dv_uniformity_json, args.site_label, args.n_years,
            args.statistic,
        )
    constrained = hydrologic_cfg is not None or dv_cfg is not None

    # --- Output directory: arm subdir then variant slug ---
    cons = "dv-l2" if args.arm == "dv_uniform" else "hydro"
    extra = {}
    if args.arm == "dv_uniform":
        extra["st"] = args.statistic
    slug = moea_slug(
        mode=args.mode,
        n_years=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        ssi=args.ssi,
        cons=cons,
        extra=extra or None,
    )
    out = stage_output_dir(STAGE, DRIVER, f"{args.arm}/{slug}")
    print(f"[04/dv_uniformity_ablation] arm={args.arm} slug={slug}")

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
            print(f"[04/dv_uniformity_ablation] anti-ideal reference from {ref}")
        else:
            print(f"[04/dv_uniformity_ablation] WARNING: --anti-ideal-reference "
                  f"{ref} not found; falling back to historical-max D*.")

    anti_ideal = compute_ssi_anti_ideal(
        hist_chars, objective_keys, headroom=cfg.anti_ideal_headroom,
        feasible_maxes=feasible_maxes,
    )
    print(f"[04/dv_uniformity_ablation] anti-ideal: {anti_ideal}")

    # --- Generator ---
    kirsch_gen = build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(
        kirsch_gen, mode=args.mode, n_years_out=args.n_years,
    )

    # --- Sanity: n_dvs vs calibration ---
    if dv_cfg is not None and generator.n_dvs != dv_cfg.n_dvs:
        print(f"[04/dv_uniformity_ablation] WARNING: generator n_dvs="
              f"{generator.n_dvs} differs from calibrated dv_cfg.n_dvs="
              f"{dv_cfg.n_dvs}.")

    # --- Config dump ---
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE,
        "driver": DRIVER,
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
    print(f"[04/dv_uniformity_ablation] arm={args.arm} Pareto: {n_pareto} solutions")

    # MM Borg: only master rank holds Pareto solutions; workers must not
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
        print(f"[04/dv_uniformity_ablation] arm={args.arm}: 0 Pareto solutions; "
              f"skipping output writes.")


if __name__ == "__main__":
    main()
