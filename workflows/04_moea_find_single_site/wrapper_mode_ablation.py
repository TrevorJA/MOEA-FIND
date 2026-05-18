"""Kirsch wrapper-mode ablation (SI-F).

Runs a single MOEA-FIND experiment under one of two Kirsch wrapper modes
so the two arms can be compared apples-to-apples:

    --wrapper-mode index    : KirschBorgWrapper with mode="index"
    --wrapper-mode residual : KirschBorgWrapper with mode="residual"

Both modes use the SAME constraint regime (DV-uniformity, Anderson-Darling
statistic, loaded from outputs/02_calibration/dv_uniformity_calibration/).
Other settings (objectives, anti-ideal, DV length, seed, NFE, algorithm)
are identical across modes; the only difference is the wrapper mode.

Compute only. Outputs to
``outputs/04_moea_find_single_site/wrapper_mode_ablation/<mode>/<slug>/``.
Per-mode diagnostic figures and the cross-mode comparison are produced by
the paired plotting drivers under ``plots/``.
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
from src.optimization.constraint_loaders import load_dv_uniformity_constraints  # noqa: E402
from src.experiment import (  # noqa: E402
    prepare_data,
    compute_historical_ssi_chars,
    compute_ssi_anti_ideal,
    extract_pareto_maxes,
    run_experiment,
)
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.optimization.constraints_dv import DVUniformityConfig, VALID_STATISTICS  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.io_paths.slugs import moea_slug  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "wrapper_mode_ablation"
VALID_WRAPPER_MODES = ("index", "residual")


def main():
    cfg = DEFAULT_EXPERIMENT
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--wrapper-mode", required=True, choices=VALID_WRAPPER_MODES)
    p.add_argument("--statistic", default="ad", choices=VALID_STATISTICS)
    p.add_argument("--nfe", type=int, default=cfg.nfe)
    p.add_argument("--n-years", type=int, default=cfg.n_years_out)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--ssi", type=int, default=cfg.ssi_timescale, choices=[1, 3, 6, 12])
    p.add_argument("--algorithm", default=cfg.algorithm,
                   choices=["borg_mm", "borg_serial"])
    p.add_argument("--dv-uniformity-json", type=Path, required=True,
                   help="DV-uniformity calibrated_dv_tolerances.json (mode-specific path)")
    p.add_argument("--site-label", default=cfg.site_label)
    p.add_argument("--anti-ideal-reference", type=Path,
                   default=cfg.anti_ideal_reference_json)
    p.add_argument("--n-islands", type=int, default=cfg.n_islands)
    p.add_argument("--checkpoint-freq", type=int, default=cfg.checkpoint_freq)
    args = p.parse_args()

    dv_cfg: DVUniformityConfig | None = load_dv_uniformity_constraints(
        args.dv_uniformity_json, args.site_label, args.n_years, args.statistic,
    )
    constrained = dv_cfg is not None

    slug = moea_slug(
        mode=args.wrapper_mode, n_years=args.n_years, nfe=args.nfe,
        seed=args.seed, ssi=args.ssi, cons="dv-l2",
        extra={"st": args.statistic},
    )
    out = stage_output_dir(STAGE, DRIVER, f"{args.wrapper_mode}/{slug}")
    print(f"[04/wrapper_mode_ablation] mode={args.wrapper_mode} slug={slug}")

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    objective_keys = cfg.objective_keys
    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(monthly_1d, args.ssi)
    feasible_maxes = None
    if args.anti_ideal_reference is not None:
        ref = Path(args.anti_ideal_reference)
        if ref.exists():
            feasible_maxes = extract_pareto_maxes(ref, objective_keys)
            print(f"  anti-ideal reference from {ref}; Pareto maxes: {feasible_maxes}")
        else:
            print(f"  WARNING: --anti-ideal-reference {ref} not found; falling back.")

    anti_ideal = compute_ssi_anti_ideal(
        hist_chars, objective_keys, headroom=cfg.anti_ideal_headroom,
        feasible_maxes=feasible_maxes,
    )
    print(f"  anti-ideal: {anti_ideal}")

    kirsch_gen = build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(kirsch_gen, mode=args.wrapper_mode,
                                  n_years_out=args.n_years)
    if dv_cfg is not None and generator.n_dvs != dv_cfg.n_dvs:
        print(f"  WARNING: generator n_dvs={generator.n_dvs} != calibrated "
              f"dv_cfg.n_dvs={dv_cfg.n_dvs}.")

    (out / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "variant": slug, "wrapper_mode": args.wrapper_mode,
        "statistic": args.statistic, "algorithm": args.algorithm,
        "nfe": args.nfe, "n_years": args.n_years, "seed": args.seed,
        "ssi": args.ssi, "constrained": constrained,
        "dv_uniformity_json": str(args.dv_uniformity_json),
        "n_islands": args.n_islands,
    }, indent=2))

    algo_kwargs = {}
    if args.algorithm in ("borg_mm", "borg_serial"):
        algo_kwargs["checkpoint_freq"] = args.checkpoint_freq
    if args.algorithm == "borg_mm":
        algo_kwargs["n_islands"] = args.n_islands

    result = run_experiment(
        monthly_2d=monthly_2d, monthly_1d=monthly_1d,
        generator=generator, generator_name=f"Kirsch ({args.wrapper_mode})",
        objective_keys=objective_keys, anti_ideal=anti_ideal,
        timescale=args.ssi, n_years_out=args.n_years,
        nfe=args.nfe, seed=args.seed, algorithm=args.algorithm,
        dv_constraint_cfg=dv_cfg, output_dir=out, **algo_kwargs,
    )
    n_pareto = result.get("n_pareto", 0)
    print(f"  Pareto: {n_pareto} solutions")

    if n_pareto > 0:
        (out / "results.json").write_text(json.dumps(result, indent=2, default=str))
        np.savez(out / "pareto.npz",
                 dvs=np.array(result["pareto_dvs"]),
                 objs=np.array(result["drought_metrics"]))


if __name__ == "__main__":
    main()
