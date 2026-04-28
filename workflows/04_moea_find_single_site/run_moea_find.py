"""Script 04 — Single-site Kirsch MOEA-FIND (manuscript Fig 5-6).

Couples Borg MOEA (or EpsNSGAII dev fallback) to the validated SynHydro
Kirsch-Nowak generator via KirschBorgWrapper. Supports both DV injection
modes (index, residual), SSI-based drought objectives, and two constraint
regimes:

    --constraint-mode dv_uniform  (production default)
        A single DV-space uniformity constraint using the Anderson-Darling
        statistic. Calibrated via bootstrap U[0,1] draws so the tolerance
        allows any configuration reachable from a uniform DV distribution.
        Chosen over the hydrologic 5-constraint set because it yields
        comparable drought-space coverage with tighter, more interpretable
        constraint geometry (see design_decisions.md §DD-13).

    --constraint-mode hydrologic  (retained for SI ablation comparison)
        Five-statistic plausibility formulation (annual mean, annual CV,
        lag-1 AC, non-drought mean, seasonal cycle) calibrated against
        historical Cannonsville flows. Used in the exp13/14 ablation.

    --constraint-mode none
        No constraints. Not used in any production or SI run.

Algorithm selection order:
    1. ``--algorithm`` CLI arg
    2. ``MOEA_FIND_ALGORITHM`` env var  (set by slurm scripts)
    3. Default: ``eps_nsga2`` (safe for local testing)

Run locally (serial EpsNSGAII, unconstrained, 20 000 NFE):
    python workflows/04_moea_find_single_site/run_moea_find.py --nfe 20000 \\
        --constraint-mode none --plot

Run with AD constraint (matches production):
    python workflows/04_moea_find_single_site/run_moea_find.py --nfe 20000 \\
        --constraint-mode dv_uniform --statistic ad

Run on HPC (MM Borg via MPI, AD constraint, 200 000 NFE):
    sbatch --export=ALL,BORG_NFE=200000 workflows/04_moea_find_single_site/slurm/run_moea_find.slurm
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

OUTPUT_SLUG = "exp04_kirsch_single_site"


_YAML_TO_CLI = {
    # YAML key → argparse dest. Two-stage parsing so --config loads YAML
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
    """Load a YAML config, returning argparse-compatible defaults.

    Unknown keys are reported and ignored so a typo in the YAML does
    not silently produce a misconfigured run.
    """
    import yaml

    raw = yaml.safe_load(Path(config_path).read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{config_path}: top-level YAML must be a mapping")
    overrides = {}
    for k, v in raw.items():
        if k in _YAML_TO_CLI:
            overrides[_YAML_TO_CLI[k]] = v
        else:
            print(f"[run_moea_find] WARN: unknown YAML key {k!r} in "
                  f"{config_path}; ignored")
    return overrides


def main():
    cfg = DEFAULT_EXPERIMENT
    # Pre-parse just --config so we can use YAML values as defaults below.
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None,
                     help="YAML preset under workflows/04_moea_find_single_site/configs/.")
    pre_args, _ = pre.parse_known_args()
    yaml_defaults = _load_yaml_config(pre_args.config) if pre_args.config else {}

    p = argparse.ArgumentParser(parents=[pre],
                                description=__doc__.splitlines()[0])
    # CLI args default to the centralized ExperimentConfig (or to YAML
    # values if --config is given). CLI flags still override either.
    p.add_argument("--nfe", type=int, default=cfg.nfe)
    p.add_argument("--n-years", type=int, default=cfg.n_years_out)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--ssi", type=int, default=cfg.ssi_timescale,
                   choices=[1, 3, 6, 12])
    p.add_argument("--mode", choices=["index", "residual"], default=cfg.dv_mode)
    p.add_argument("--algorithm", default=cfg.algorithm,
                   choices=["borg_mm", "borg_serial", "eps_nsga2"],
                   help="MOEA backend. Overridden by MOEA_FIND_ALGORITHM env var.")
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
                        "drives D* placement for NON-cyclic objectives. "
                        "Cyclic metrics still use 12*headroom. If unset, "
                        "D* falls back to historical_max*headroom.")
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
    p.add_argument("--plot", action="store_true")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    if yaml_defaults:
        p.set_defaults(**yaml_defaults)
    args = p.parse_args()

    # --- Constraints (load early — needed for variant slug) ---
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
    out = args.output_dir / slug
    out.mkdir(parents=True, exist_ok=True)
    print(f"[04] variant: {slug}")

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
    # Augment historical chars with the trace-level extras (Q10 flow,
    # time-in-drought) so metric extractors that read those keys work.
    from src.objectives import compute_ssi_drought_characteristics
    hist_chars = compute_ssi_drought_characteristics(
        ssi_hist, monthly_flows=monthly_1d
    )
    feasible_maxes = None
    if args.anti_ideal_reference is not None:
        ref = Path(args.anti_ideal_reference)
        if ref.exists():
            feasible_maxes = extract_pareto_maxes(ref, metric_set)
            print(f"[04] anti-ideal reference from {ref}")
            print(f"     Pareto maxes: {feasible_maxes}")
        else:
            print(f"[04] WARNING: --anti-ideal-reference {ref} not found; "
                  f"falling back to historical-max D*.")

    anti_ideal = compute_ssi_anti_ideal(
        hist_chars,
        metric_set,
        headroom=cfg.anti_ideal_headroom,
        feasible_maxes=feasible_maxes,
    )
    print(f"[04] historical SSI-{args.ssi}: n={hist_chars['n_events']}")
    for m in metric_set:
        print(f"     {m.label}: {m.extract(hist_chars):.4g} ({m.units})")
    print(f"[04] metric set: {args.metric_set} → {objective_keys}")
    print(f"[04] anti-ideal: {anti_ideal}")

    # --- Generator ---
    kirsch_gen = build_kirsch_generator(monthly_2d)
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

    print(f"[04] Pareto: {result.get('n_pareto', 0)} solutions")

    # Under MM Borg MPI, only the master rank has Pareto solutions.
    # Workers report 0 and must NOT overwrite the master's output files.
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
    else:
        print(f"[04] WARNING: 0 Pareto solutions with {args.nfe} NFE "
              f"and {generator.n_dvs} DVs. Skipping output writes "
              f"(worker rank under MM Borg MPI, or no solutions found).")

    # --- Plots (saved to variant directory, not global figures/) ---
    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from src.plotting.drought_space import (
                plot_scatter_with_marginals,
                plot_drought_space_3d,
            )
            from src.plotting.trace_diagnostics import (
                plot_autocorrelation_comparison,
                plot_flow_duration_curve,
                plot_hydrology_panels,
                plot_seasonal_cycle_comparison,
            )

            fig_dir = out / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)

            pareto_chars = result.get("pareto_chars", [])
            if pareto_chars:
                dm = np.array(result["drought_metrics"])
                axis_labels = tuple(
                    f"{m.label} ({m.units})" for m in metric_set
                )
                hist_point = tuple(
                    float(m.extract(hist_chars)) for m in metric_set[:2]
                )
                # Compute per-block historical drought characteristics
                # using the SAME prefitted SSI calculator the MOEA used,
                # so the historical cloud is directly comparable to the
                # Pareto archive in drought-characteristic space.
                from src.historical_blocks import compute_historical_block_chars
                hist_block_chars = compute_historical_block_chars(
                    monthly_1d,
                    T_years=args.n_years,
                    ssi_calc=ssi_calc,
                    objective_keys=objective_keys,
                    stride=1,
                )
                print(f"[04] historical block drought-chars: "
                      f"{hist_block_chars.shape[0]} blocks, "
                      f"{objective_keys[0]} ["
                      f"{hist_block_chars[:, 0].min():.2f}, "
                      f"{hist_block_chars[:, 0].max():.2f}]")
                fig_a = plot_scatter_with_marginals(
                    dm[:, :2],
                    title=f"Kirsch ({args.mode}) Pareto vs historical blocks",
                    historical_point=hist_point,
                    anti_ideal=anti_ideal[:2],
                    historical_cloud=hist_block_chars[:, :2],
                    objective_labels=axis_labels[:2],
                )
                fig_a.savefig(fig_dir / "fig05_drought_space.pdf", dpi=300)
                plt.close(fig_a)
                print(f"[04] wrote {fig_dir / 'fig05_drought_space.pdf'}")

                # Stash the block chars for downstream / 3D plotting.
                np.savez(out / "historical_block_chars.npz",
                         chars=hist_block_chars,
                         objective_keys=np.array(list(objective_keys)))

                # 3D drought-space scatter (if three objectives were used).
                if dm.shape[1] >= 3 and len(metric_set) >= 3:
                    hist_point_3d = np.array([
                        float(m.extract(hist_chars)) for m in metric_set[:3]
                    ])
                    fig_3d = plot_drought_space_3d(
                        dm[:, :3],
                        anti_ideal=anti_ideal[:3],
                        objective_labels=axis_labels[:3],
                        historical_point=hist_point_3d,
                        historical_cloud=hist_block_chars[:, :3],
                        title=(
                            f"Kirsch ({args.mode}) Pareto — "
                            f"NFE={args.nfe}, n={dm.shape[0]}"
                        ),
                    )
                    fig_3d.savefig(fig_dir / "fig05_drought_space_3d.pdf", dpi=200,
                                    bbox_inches="tight")
                    fig_3d.savefig(fig_dir / "fig05_drought_space_3d.png", dpi=200,
                                    bbox_inches="tight")
                    plt.close(fig_3d)
                    print(f"[04] wrote {fig_dir / 'fig05_drought_space_3d.pdf'}")

            pareto_traces_1d = result.get("pareto_traces_1d", [])
            pareto_traces_2d = result.get("pareto_traces_2d", [])
            if pareto_traces_1d:
                traces_1d = [np.array(t) for t in pareto_traces_1d]
                traces_2d = [np.array(t) for t in pareto_traces_2d]

                # Build historical block ensemble at the same length as
                # the synthetic traces. This is the fair comparator for
                # FDC/ACF/seasonal diagnostics.
                from src.historical_blocks import (
                    resample_historical_blocks,
                    resample_historical_blocks_2d,
                )
                hist_blocks_1d = resample_historical_blocks(
                    monthly_1d, T_years=args.n_years, stride=1,
                )
                hist_blocks_2d = resample_historical_blocks_2d(
                    monthly_2d, T_years=args.n_years, stride=1,
                )
                print(f"[04] historical blocks: "
                      f"{len(hist_blocks_1d)} overlapping {args.n_years}-yr blocks")

                fig_acf, _ = plot_autocorrelation_comparison(
                    traces_1d, hist_blocks_1d,
                )
                fig_acf.savefig(fig_dir / "fig06a_acf.pdf", dpi=300)
                plt.close(fig_acf)

                fig_fdc, _ = plot_flow_duration_curve(
                    traces_1d, hist_blocks_1d,
                )
                fig_fdc.savefig(fig_dir / "fig06b_fdc.pdf", dpi=300)
                plt.close(fig_fdc)

                fig_sc, _ = plot_seasonal_cycle_comparison(
                    traces_2d, hist_blocks_2d,
                )
                fig_sc.savefig(fig_dir / "fig06c_seasonal.pdf", dpi=300)
                plt.close(fig_sc)

                fig_hy, _ = plot_hydrology_panels(
                    traces_1d, hist_blocks_1d, traces_2d, hist_blocks_2d,
                )
                fig_hy.savefig(fig_dir / "fig05_hydrology.pdf", dpi=300,
                               bbox_inches="tight")
                plt.close(fig_hy)

                print(f"[04] wrote trace diagnostics to {fig_dir}")

        except ImportError as exc:
            print(f"[04] skipping plots (import error: {exc})")


if __name__ == "__main__":
    main()
