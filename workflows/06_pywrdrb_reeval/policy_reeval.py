"""Stage 06 / policy_reeval -- DRB policy re-evaluation (manuscript Fig 9 data).

Consumes a MOEA-FIND Pareto archive from stage 04 (single-site) or stage 05
(multi-site -- now folded into this driver via
``replay_pareto_to_multisite_monthly``). Replays Pareto DVs through a
multi-site Kirsch-Nowak-KDE pipeline, converts to Pywr-DRB FlowEnsemble HDF5
format, runs Pywr-DRB simulations with constant_max demand, extracts FFMP
drought levels, and produces a per-realization metric bank.

Spatial correlation across the multi-site ensemble is preserved by
construction: ``replay_pareto_to_multisite_monthly`` replays each Pareto
DV vector through *shared* monthly Kirsch indexes, so all sites see the
same monthly resampling pattern. No separate stage 05 driver is needed.

Four pipeline stages (run independently via ``--stage``):
    generate  -- Replay DVs -> multi-site daily flows -> HDF5
    prep      -- Run pywrdrb predicted-inflow preprocessor
    simulate  -- Run Pywr-DRB model for each Pareto realization
    classify  -- Extract Level 6, build satisficing table, write metric bank

Outputs under ``outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/``:
    - config.json
    - pywrdrb_inputs/{gage_flow_mgd.hdf5, catchment_inflow_mgd.hdf5,
                      predicted_inflows_mgd.hdf5}
    - simulations/pywrdrb_output.hdf5
    - results/{metric_bank.parquet, satisficing_table.parquet,
               drought_levels.npz}

Plotting / scenario-discovery figures live in stage 07.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.multisite_data import (  # noqa: E402
    load_pywrdrb_gage_flow,
    load_pywrdrb_catchment_inflow,
    get_kirsch_sites,
    get_kde_regression_sites,
    get_kde_pairs,
    fit_multisite_generators,
)
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402
from src.pywrdrb_bridge import (  # noqa: E402
    replay_pareto_to_multisite_monthly,
    disaggregate_monthly_to_daily,
    fit_kde_models,
    generate_kde_downstream_nodes,
    compute_marginal_catchment_inflows,
    write_flowensemble_hdf5,
    register_flow_type,
    prep_predicted_inflows,
    run_pywrdrb_batch,
)
from src.scenario_discovery import (  # noqa: E402
    extract_drought_levels,
    build_satisficing_table,
    save_results,
)
from src.satisficing_metrics import (  # noqa: E402
    compute_metric_bank,
    write_metric_bank,
)

STAGE = "06_pywrdrb_reeval"
DRIVER = "policy_reeval"

VALID_STAGES = ("all", "generate", "prep", "simulate", "classify")


def main():
    p = argparse.ArgumentParser(
        description="DRB policy re-evaluation: MOEA-FIND Pareto -> Pywr-DRB"
    )
    p.add_argument(
        "--pareto-results", type=Path, required=True,
        help="Path to results.json from stage 04 (run_moea_find/<slug>/results.json).",
    )
    p.add_argument(
        "--stage", default="all", choices=VALID_STAGES,
        help="Pipeline stage to run (default: all).",
    )
    p.add_argument("--n-years", type=int, default=20,
                   help="Synthetic trace length (must match stage 04).")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for stage 04 results.")
    p.add_argument("--mode", default="residual",
                   choices=["index", "residual"],
                   help="DV injection mode (must match stage 04).")
    p.add_argument("--start-date", default="2030-01-01",
                   help="Pywr-DRB simulation start date.")
    p.add_argument("--end-date", default=None,
                   help="Pywr-DRB simulation end date (default: computed "
                        "from n-years and start-date).")
    p.add_argument("--baseline-dataset",
                   default="pub_nhmv10_BC_withObsScaled",
                   help="Pywr-DRB flow dataset for multi-site fitting.")
    p.add_argument("--demand-source", default="constant_max",
                   help="NYC/NJ demand source for Pywr-DRB.")
    p.add_argument("--batch-size", type=int, default=10,
                   help="Realizations per Pywr-DRB simulation batch.")
    p.add_argument("--disagg-seed", type=int, default=42,
                   help="Seed for Nowak disaggregation.")
    p.add_argument("--use-mpi", action="store_true",
                   help="Parallelize Stages 2-3 across MPI ranks.")
    p.add_argument("--realization-subset-file", type=Path, default=None,
                   help="Optional JSON listing a subset of realization indices "
                        "(integer list) to process. Used for dev-ensemble smoke "
                        "tests -- the production run omits this flag and "
                        "processes all Pareto solutions.")
    args = p.parse_args()

    # --- Load Pareto results from stage 04 ---
    print(f"[06/policy_reeval] Loading Pareto results from {args.pareto_results}")
    results = json.loads(args.pareto_results.read_text())
    pareto_dvs = np.array(results["pareto_dvs"])
    pareto_chars = results.get("pareto_chars", [])
    drought_metrics = np.array(results.get("drought_metrics", []))
    objective_keys = tuple(results.get("objective_keys",
                                       ["mean_duration", "mean_avg_severity"]))
    n_pareto = results.get("n_pareto", len(pareto_dvs))

    print(f"[06/policy_reeval] {n_pareto} Pareto solutions, "
          f"{len(pareto_dvs[0])} DVs, objectives={objective_keys}")

    # Infer n_years from DVs if not specified
    n_years = args.n_years
    if args.mode == "residual" and len(pareto_dvs[0]) != n_years * 12:
        inferred = len(pareto_dvs[0]) // 12
        print(f"[06/policy_reeval] WARNING: --n-years={n_years} but DVs imply "
              f"T={inferred}. Using T={inferred}.")
        n_years = inferred

    # Compute end date from start date and n_years
    if args.end_date is None:
        from dateutil.relativedelta import relativedelta
        from datetime import datetime
        start_dt = datetime.strptime(args.start_date, "%Y-%m-%d")
        end_dt = start_dt + relativedelta(years=n_years) - relativedelta(days=1)
        end_date = end_dt.strftime("%Y-%m-%d")
    else:
        end_date = args.end_date

    # --- Output directory keyed by SOURCE slug (parent dir of results.json) ---
    # The source slug embeds every knob set by stage 04 (mode, T, NFE, seed,
    # constraint flag, constraint mode, statistic). Reusing it here keeps the
    # 06 outputs paired 1:1 with their 04 archive on disk.
    src_slug = args.pareto_results.parent.name
    out = stage_output_dir(STAGE, DRIVER, src_slug)

    pywrdrb_inputs = out / "pywrdrb_inputs"
    sim_dir = out / "simulations"
    results_dir = out / "results"

    flow_type = f"moea_find_{src_slug}"

    # Apply --realization-subset-file if present. The subset file holds
    # *global* Pareto indices; we slice the archive down to the selected
    # rows, then re-number to *local* indices 0..N-1 for every Pywr-DRB-
    # facing data structure. ``replay_pareto_to_multisite_monthly`` writes
    # the catchment-inflow HDF5 with keys 0..N-1 (enumerate order), so
    # ``realization_ids`` consumed by the preprocessors must match that
    # local numbering. The global indices are preserved in config.json
    # for traceability back to stage 04's Pareto ordering.
    global_pareto_indices: List[int]
    if args.realization_subset_file is not None:
        subset_path = Path(args.realization_subset_file)
        subset_idx = json.loads(subset_path.read_text())
        subset_idx = [int(i) for i in subset_idx]
        if any(i < 0 or i >= n_pareto for i in subset_idx):
            raise SystemExit(
                f"[06/policy_reeval] realization subset contains out-of-range "
                f"indices (n_pareto={n_pareto}): {subset_idx}"
            )
        pareto_dvs = pareto_dvs[subset_idx]
        pareto_chars = [pareto_chars[i] for i in subset_idx] if pareto_chars else []
        drought_metrics = drought_metrics[subset_idx] if drought_metrics.size else drought_metrics
        global_pareto_indices = subset_idx
        n_pareto = len(subset_idx)
        print(f"[06/policy_reeval] realization subset: {n_pareto} members "
              f"from {subset_path} (global Pareto indices {subset_idx})")
    else:
        global_pareto_indices = list(range(n_pareto))

    # Local indices 0..N-1 -- what pywrdrb preprocessors + simulators see.
    realization_ids = [str(i) for i in range(n_pareto)]

    # --- Config dump ---
    config = {
        "stage": STAGE,
        "driver": DRIVER,
        "pareto_source": str(args.pareto_results),
        "src_slug": src_slug,
        "flow_type": flow_type,
        "n_pareto": n_pareto,
        "n_years": n_years,
        "mode": args.mode,
        "seed": args.seed,
        "start_date": args.start_date,
        "end_date": end_date,
        "baseline_dataset": args.baseline_dataset,
        "demand_source": args.demand_source,
        "batch_size": args.batch_size,
        "objective_keys": list(objective_keys),
        "constrained": results.get("constraint_config") is not None,
        "global_pareto_indices": global_pareto_indices,
    }
    (out / "config.json").write_text(json.dumps(config, indent=2))
    print(f"[06/policy_reeval] src_slug: {src_slug}")
    print(f"[06/policy_reeval] output: {out}")
    print(f"[06/policy_reeval] simulation period: {args.start_date} to {end_date}")

    # ===============================================================
    # Stage 1: GENERATE -- multi-site daily flows from Pareto DVs
    # ===============================================================
    if args.stage in ("all", "generate"):
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"[06/policy_reeval] STAGE 1: GENERATE")
        print(f"{'='*60}")

        # Load multi-site historical data
        Q_gage = load_pywrdrb_gage_flow(args.baseline_dataset)
        Q_inflow = load_pywrdrb_catchment_inflow(args.baseline_dataset)

        kirsch_sites = get_kirsch_sites(Q_gage)
        kde_sites = get_kde_regression_sites(Q_gage)
        kde_pairs = get_kde_pairs(kirsch_sites, kde_sites)

        print(f"[06/policy_reeval] Kirsch sites: {len(kirsch_sites)}, "
              f"KDE sites: {len(kde_sites)}, "
              f"KDE pairs: {len(kde_pairs)}")

        # Fit multi-site generators
        kirsch_gen, nowak_disagg = fit_multisite_generators(
            Q_gage, kirsch_sites
        )

        # Wrap with KirschBorgWrapper
        wrapper = KirschBorgWrapper(
            kirsch_gen, mode=args.mode, n_years_out=n_years
        )
        print(f"[06/policy_reeval] Wrapper: {wrapper.n_dvs} DVs, "
              f"{wrapper.n_sites} sites, mode={args.mode}")

        # Replay DVs -> multi-site monthly (shared monthly Kirsch indexes
        # across sites preserve spatial correlation by construction).
        monthly = replay_pareto_to_multisite_monthly(
            pareto_dvs, wrapper, kirsch_sites,
            start_date=args.start_date,
        )

        # Monthly -> daily (Nowak)
        daily = disaggregate_monthly_to_daily(
            monthly, nowak_disagg, seed=args.disagg_seed,
        )

        # KDE downstream nodes
        kdes = fit_kde_models(Q_inflow, kde_pairs)
        daily_full = generate_kde_downstream_nodes(
            daily, kdes, kde_pairs, seed=args.disagg_seed,
        )

        # Save gage flow HDF5 (needed for diversion preprocessing if used)
        pywrdrb_inputs.mkdir(parents=True, exist_ok=True)
        write_flowensemble_hdf5(
            daily_full, pywrdrb_inputs / "gage_flow_mgd.hdf5",
        )

        # Gage flows -> marginal catchment inflows
        catchment_inflows = compute_marginal_catchment_inflows(daily_full)

        # Write catchment inflow HDF5 (FlowEnsemble format)
        write_flowensemble_hdf5(
            catchment_inflows, pywrdrb_inputs / "catchment_inflow_mgd.hdf5",
        )

        elapsed = time.time() - t0
        print(f"[06/policy_reeval] Stage 1 complete in {elapsed:.1f}s")

    # ===============================================================
    # Stage 2: PREP -- pywrdrb preprocessors
    # ===============================================================
    if args.stage in ("all", "prep"):
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"[06/policy_reeval] STAGE 2: PREP")
        print(f"{'='*60}")

        register_flow_type(flow_type, pywrdrb_inputs)
        prep_predicted_inflows(
            flow_type, pywrdrb_inputs, realization_ids,
            use_mpi=args.use_mpi,
            demand_source=args.demand_source,
        )

        elapsed = time.time() - t0
        print(f"[06/policy_reeval] Stage 2 complete in {elapsed:.1f}s")

    # ===============================================================
    # Stage 3: SIMULATE -- Pywr-DRB model runs
    # ===============================================================
    if args.stage in ("all", "simulate"):
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"[06/policy_reeval] STAGE 3: SIMULATE")
        print(f"{'='*60}")

        # Re-register in case running stage independently
        register_flow_type(flow_type, pywrdrb_inputs)

        run_pywrdrb_batch(
            flow_type=flow_type,
            realization_ids=realization_ids,
            start_date=args.start_date,
            end_date=end_date,
            output_dir=sim_dir,
            model_dir=out / "models",
            batch_size=args.batch_size,
            demand_source=args.demand_source,
            use_mpi=args.use_mpi,
        )

        elapsed = time.time() - t0
        print(f"[06/policy_reeval] Stage 3 complete in {elapsed:.1f}s")

    # ===============================================================
    # Stage 4: CLASSIFY -- metric bank + legacy FFMP table
    # ===============================================================
    if args.stage in ("all", "classify"):
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"[06/policy_reeval] STAGE 4: CLASSIFY")
        print(f"{'='*60}")

        results_dir.mkdir(parents=True, exist_ok=True)

        # Locate the combined Pywr-DRB output (produced by Stage 3)
        combined_hdf5 = sim_dir / "pywrdrb_output.hdf5"
        if not combined_hdf5.exists():
            # Tolerate running Stage 4 separately by scanning for
            # rank/batch leftovers as a fallback.
            candidates = sorted(sim_dir.glob("*.hdf5"))
            if candidates:
                combined_hdf5 = candidates[0]
                print(f"[06/policy_reeval] combined output not found; falling "
                      f"back to {combined_hdf5}")
            else:
                raise FileNotFoundError(
                    f"No Pywr-DRB output HDF5 found under {sim_dir}"
                )

        # Metric bank -- always computed; cheap compared to simulation
        print(f"[06/policy_reeval] computing metric bank from {combined_hdf5}")
        bank = compute_metric_bank(combined_hdf5, realization_ids)
        bank_path = results_dir / "metric_bank.parquet"
        bank_written_to = write_metric_bank(bank, bank_path)
        print(f"[06/policy_reeval] wrote metric bank: {bank_written_to} "
              f"({len(bank)} realizations, {bank.shape[1]} metrics)")

        # Legacy FFMP-Level-6 classification table. We still compute and
        # persist this because it's the simplest drought-level summary
        # that downstream plotting (stage 07) and any analyst needs
        # before invoking a classifier. No figures are generated here;
        # all scenario-discovery plotting lives in stage 07.
        drought_levels = extract_drought_levels(sim_dir, realization_ids)
        df = build_satisficing_table(
            drought_levels, pareto_chars, drought_metrics, objective_keys,
        )
        save_results(df, drought_levels, results_dir)

        elapsed = time.time() - t0
        print(f"[06/policy_reeval] Stage 4 complete in {elapsed:.1f}s")

    print(f"\n[06/policy_reeval] Pipeline complete. Results in {out}")


if __name__ == "__main__":
    main()
