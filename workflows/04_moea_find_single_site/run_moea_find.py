"""Stage 04 / run_moea_find -- single-site Kirsch MOEA-FIND (Figs 5-6).

Couples MM Borg MOEA to the validated
SynHydro Kirsch-Nowak generator via :class:`KirschBorgWrapper`. Writes
only numerical artifacts; figures are produced by the paired plotting
driver ``workflows/04_moea_find_single_site/plots/run_moea_find.py``.

Outputs under ``outputs/04_moea_find_single_site/run_moea_find/<slug>/``:
    config.json
    results.json
    pareto.npz
    historical_block_chars.npz   (per-block historical drought-chars)
    historical_blocks.npz        (resampled historical 1d/2d blocks)

Constraint regimes:
    --constraint-mode dv_uniform  (production default)
        Single DV-space uniformity constraint, Anderson-Darling statistic.
        Calibration JSON path: outputs/02_calibration/dv_uniformity_calibration/<mode>/calibrated_dv_tolerances.json
    --constraint-mode hydrologic  (SI ablation comparison)
        Five-statistic plausibility formulation.
    --constraint-mode none        (development only)
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
from src.drought_metrics import (  # noqa: E402
    PRESETS,
    metric_labels,
    metric_names,
    resolve_metric_set,
)
from src.paths import stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "run_moea_find"


_YAML_TO_CLI = {
    # YAML key -> argparse dest. Two-stage parsing so --config loads YAML
    # values into the parser defaults, then a second parse lets CLI flags
    # override.
    "nfe": "nfe",
    "n_years": "n_years",
    "seed": "seed",
    "ssi": "ssi",
    "mode": "mode",
    "algorithm": "algorithm",
    "constraint_mode": "constraint_mode",
    "statistic": "statistic",
    "dv_uniformity_json": "dv_uniformity_json",
    "constraints_json": "constraints_json",
    "anti_ideal_reference": "anti_ideal_reference",
    "site_label": "site_label",
    "metric_set": "metric_set",
    "n_islands": "n_islands",
    "checkpoint_freq": "checkpoint_freq",
    "old_checkpoint": "old_checkpoint",
}


def _load_yaml_config(config_path: Path) -> dict:
    """Load a YAML config, returning argparse-compatible defaults."""
    import yaml

    raw = yaml.safe_load(Path(config_path).read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{config_path}: top-level YAML must be a mapping")
    overrides = {}
    for k, v in raw.items():
        if k in _YAML_TO_CLI:
            overrides[_YAML_TO_CLI[k]] = v
        else:
            print(f"[04/run_moea_find] WARN: unknown YAML key {k!r} in "
                  f"{config_path}; ignored")
    return overrides


def main():
    cfg = DEFAULT_EXPERIMENT
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None,
                     help="YAML preset under workflows/04_moea_find_single_site/configs/.")
    pre_args, _ = pre.parse_known_args()
    yaml_defaults = _load_yaml_config(pre_args.config) if pre_args.config else {}

    p = argparse.ArgumentParser(parents=[pre],
                                description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=cfg.nfe)
    p.add_argument("--n-years", type=int, default=cfg.n_years_out)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--ssi", type=int, default=cfg.ssi_timescale,
                   choices=[1, 3, 6, 12])
    p.add_argument("--mode", choices=["index", "residual"], default=cfg.dv_mode)
    p.add_argument("--algorithm", default=cfg.algorithm,
                   choices=["borg_mm", "borg_serial"],
                   help="Borg backend. Overridden by MOEA_FIND_ALGORITHM env var.")
    p.add_argument("--constraint-mode",
                   choices=["dv_uniform", "hydrologic", "none"],
                   default=cfg.constraint_mode,
                   help="Which constraint regime to use. Production default: dv_uniform.")
    p.add_argument("--statistic", default=cfg.dv_uniformity_statistic,
                   choices=VALID_STATISTICS,
                   help="DV-uniformity statistic (only used when --constraint-mode=dv_uniform).")
    p.add_argument("--dv-uniformity-json", type=Path, default=cfg.dv_uniformity_json,
                   help="Path to calibrated_dv_tolerances.json (dv_uniform mode).")
    p.add_argument("--constraints-json", type=Path, default=cfg.constraints_json,
                   help="Path to calibrated_tolerances.json (hydrologic mode).")
    p.add_argument("--anti-ideal-reference", type=Path,
                   default=cfg.anti_ideal_reference_json,
                   help="Path to a prior results.json whose Pareto max "
                        "drives D* placement for non-cyclic objectives.")
    p.add_argument("--site-label", default=cfg.site_label)
    p.add_argument(
        "--metric-set", default=cfg.metric_set,
        help=(
            "Drought metric set: a preset name "
            f"({', '.join(sorted(PRESETS.keys()))}) or a single metric "
            "name from the registry."
        ),
    )
    p.add_argument("--n-islands", type=int, default=cfg.n_islands,
                   help="Number of islands for MM Borg.")
    p.add_argument("--checkpoint-freq", type=int, default=cfg.checkpoint_freq)
    p.add_argument("--old-checkpoint", type=Path, default=None,
                   help="Path to a checkpoint file to resume from.")
    if yaml_defaults:
        p.set_defaults(**yaml_defaults)
    args = p.parse_args()

    # --- Constraints (load early -- needed for variant slug) ---
    constraint_cfg = None
    dv_constraint_cfg = None
    if args.constraint_mode == "dv_uniform":
        dv_constraint_cfg = load_dv_uniformity_constraints(
            args.dv_uniformity_json, args.site_label, args.n_years, args.statistic,
        )
        constrained = dv_constraint_cfg is not None
    elif args.constraint_mode == "hydrologic":
        constraint_cfg = load_hydrologic_constraints(
            args.constraints_json, args.site_label, args.n_years,
        )
        constrained = constraint_cfg is not None
    else:
        constrained = False

    # --- Output directory keyed by variant slug ---
    slug_extra: dict = {}
    if args.constraint_mode == "dv_uniform":
        slug_extra["cm"] = "dv_uniform"
        slug_extra["st"] = args.statistic
    elif args.constraint_mode == "hydrologic":
        slug_extra["cm"] = "hydrologic"
    slug = make_variant_slug(
        mode=args.mode,
        n_years=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        constrained=constrained,
        extra=slug_extra if slug_extra else None,
    )
    out = stage_output_dir(STAGE, DRIVER, slug)
    print(f"[04/run_moea_find] variant: {slug}")

    # --- Data ---
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    # --- SSI characterisation ---
    metric_set = resolve_metric_set(args.metric_set)
    objective_keys = metric_names(metric_set)
    ssi_hist, ssi_calc, hist_chars = compute_historical_ssi_chars(
        monthly_1d, args.ssi
    )
    from src.objectives import compute_ssi_drought_characteristics
    hist_chars = compute_ssi_drought_characteristics(
        ssi_hist, monthly_flows=monthly_1d
    )
    feasible_maxes = None
    if args.anti_ideal_reference is not None:
        ref = Path(args.anti_ideal_reference)
        if ref.exists():
            feasible_maxes = extract_pareto_maxes(ref, metric_set)
            print(f"[04/run_moea_find] anti-ideal reference from {ref}")
            print(f"     Pareto maxes: {feasible_maxes}")
        else:
            print(f"[04/run_moea_find] WARNING: --anti-ideal-reference {ref} "
                  f"not found; falling back to historical-max D*.")

    anti_ideal = compute_ssi_anti_ideal(
        hist_chars,
        metric_set,
        headroom=cfg.anti_ideal_headroom,
        feasible_maxes=feasible_maxes,
    )
    print(f"[04/run_moea_find] historical SSI-{args.ssi}: "
          f"n={hist_chars['n_events']}")
    for m in metric_set:
        print(f"     {m.label}: {m.extract(hist_chars):.4g} ({m.units})")
    print(f"[04/run_moea_find] metric set: {args.metric_set} -> {objective_keys}")
    print(f"[04/run_moea_find] anti-ideal: {anti_ideal}")

    # --- Generator ---
    kirsch_gen = build_kirsch_generator(monthly_2d)
    generator = KirschBorgWrapper(
        kirsch_gen, mode=args.mode, n_years_out=args.n_years,
    )

    # --- Config dump ---
    (out / "config.json").write_text(json.dumps({
        "stage": STAGE,
        "driver": DRIVER,
        "variant": slug,
        "algorithm": args.algorithm,
        "nfe": args.nfe,
        "n_years": args.n_years,
        "seed": args.seed,
        "ssi": args.ssi,
        "mode": args.mode,
        "constrained": constrained,
        "constraint_mode": args.constraint_mode,
        "statistic": args.statistic if args.constraint_mode == "dv_uniform" else None,
        "constraints_json": str(args.constraints_json) if args.constraints_json else None,
        "dv_uniformity_json": str(args.dv_uniformity_json) if args.dv_uniformity_json else None,
        "n_islands": args.n_islands,
        "metric_set": args.metric_set,
        "objective_keys": list(objective_keys),
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
        objective_keys=metric_set,
        anti_ideal=anti_ideal,
        timescale=args.ssi,
        n_years_out=args.n_years,
        nfe=args.nfe,
        seed=args.seed,
        algorithm=args.algorithm,
        constraint_cfg=constraint_cfg,
        dv_constraint_cfg=dv_constraint_cfg,
        output_dir=out,
        **algo_kwargs,
    )

    print(f"[04/run_moea_find] Pareto: {result.get('n_pareto', 0)} solutions")

    # Under MM Borg MPI only the master rank holds Pareto solutions; workers
    # report 0 and must NOT overwrite the master's output files.
    if result.get("n_pareto", 0) > 0:
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

        # Cache the historical block characteristics + raw historical
        # blocks alongside the Pareto archive so the plotting driver can
        # render Figs 5/6 without re-doing the historical-block compute.
        from src.historical_blocks import (
            compute_historical_block_chars,
            resample_historical_blocks,
            resample_historical_blocks_2d,
        )
        hist_block_chars = compute_historical_block_chars(
            monthly_1d, T_years=args.n_years, ssi_calc=ssi_calc,
            objective_keys=objective_keys, stride=1,
        )
        np.savez(
            out / "historical_block_chars.npz",
            chars=hist_block_chars,
            objective_keys=np.array(list(objective_keys)),
        )
        hist_blocks_1d = resample_historical_blocks(
            monthly_1d, T_years=args.n_years, stride=1,
        )
        hist_blocks_2d = resample_historical_blocks_2d(
            monthly_2d, T_years=args.n_years, stride=1,
        )
        np.savez(
            out / "historical_blocks.npz",
            blocks_1d=np.array(hist_blocks_1d),
            blocks_2d=np.array(hist_blocks_2d),
            anti_ideal=anti_ideal,
        )
        print(f"     wrote {out / 'historical_block_chars.npz'} and "
              f"{out / 'historical_blocks.npz'}")
    else:
        print(f"[04/run_moea_find] WARNING: 0 Pareto solutions with "
              f"{args.nfe} NFE and {generator.n_dvs} DVs. Skipping output "
              f"writes (worker rank under MM Borg MPI, or no solutions found).")


if __name__ == "__main__":
    main()
