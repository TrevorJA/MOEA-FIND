"""Stage 06 (windowed) -- per-scenario Pywr-DRB SIMULATION (consumes prep).

Slim simulation-only driver after the 2026-05-20 prep/sim split. Reads
the per-scenario inputs and window metadata that ``prepare_windowed.py``
wrote, then runs Pywr-DRB over each scenario's first-drought window.
Re-running the sim (e.g. to change demand source or pick up a Pywr-DRB
bug fix) does NOT redo the heavy Kirsch / Nowak / KDE chain.

Pipeline:
    1) prepare_windowed.py        -- DV -> daily multisite -> HDF5
    2) policy_reeval_windowed.py  -- HDF5 + window -> pywrdrb sim   <-- this
    3) aggregate_windowed.py      -- per-scenario -> metric bank

Outputs (per scenario ``s`` under the COARSE src slug):
  outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/
    scenarios/<s>/
      simulations/batch_0.hdf5   (pywrdrb OutputRecorder; this script)
      status_sim.json            (ok/skip/fail)

A separate aggregation pass (``aggregate_windowed.py``) concatenates the
per-scenario hdf5 + window metadata into ``results/metric_bank.parquet``
for stages 07/08/09.

Design: EMBARRASSINGLY PARALLEL, zero shared state / MPI collectives.
Each (array task, MPI rank) pair owns a disjoint scenario subset.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.pywrdrb.bridge import (  # noqa: E402
    register_flow_type,
    _get_mpi_context,
)
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "06_pywrdrb_reeval"
DRIVER = "policy_reeval"


def _worker_scenarios(n_total: int, array_id: int, array_count: int,
                      rank: int, size: int) -> list[int]:
    """Disjoint scenario subset for (array task, MPI rank).

    Same arithmetic as prepare_windowed so the prep and sim arrays
    target identical scenario assignments at matching (array_id, rank).
    """
    task_chunk = np.array_split(np.arange(n_total), array_count)[array_id]
    return [int(i) for i in np.array_split(task_chunk, size)[rank]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pareto-results", type=Path, required=True,
                    help="Stage-04 results.json (used only for scenario "
                         "count and src_slug; DV/window data come from "
                         "the per-scenario prep outputs).")
    ap.add_argument("--site-label", default="cannonsville")
    ap.add_argument("--demand-source", default="constant_max",
                    help="Demand label for the Pywr-DRB ModelBuilder.")
    ap.add_argument("--array-id", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--array-count", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1)))
    ap.add_argument("--use-mpi", action="store_true")
    ap.add_argument("--overwrite", action="store_true",
                    help="Re-run scenarios whose status_sim.json reports "
                         "ok; default behavior skips them.")
    args = ap.parse_args()

    comm = None
    rank, size = 0, 1
    if args.use_mpi:
        comm, rank, size = _get_mpi_context()

    results = json.loads(Path(args.pareto_results).read_text())
    n_total = len(np.asarray(results["pareto_dvs"], dtype=float))
    src_slug = args.pareto_results.parent.name

    out = stage_output_dir(STAGE, DRIVER, src_slug)
    scen_root = out / "scenarios"

    my_ids = _worker_scenarios(n_total, args.array_id, args.array_count,
                               rank, size)
    print(f"[06sim] array={args.array_id}/{args.array_count} "
          f"rank={rank}/{size} -> {len(my_ids)} scenarios of {n_total} "
          f"(src={src_slug})", flush=True)
    if not my_ids:
        return

    n_ok, n_skip_done, n_skip_noprep, n_fail = 0, 0, 0, 0
    for sid in my_ids:
        t0 = time.time()
        sdir = scen_root / str(sid)
        sin = sdir / "pywrdrb_inputs"
        prep_status = sdir / "status_prep.json"
        sim_status = sdir / "status_sim.json"
        window_path = sdir / "window.json"

        # Skip if a successful sim already exists (idempotent re-runs).
        if not args.overwrite and sim_status.exists():
            try:
                if json.loads(sim_status.read_text()).get("ok"):
                    n_skip_done += 1
                    continue
            except Exception:  # noqa: BLE001
                pass

        # Require a successful prep (status_prep.ok=True) and the
        # window metadata. Missing prep is a hard skip; rerun
        # prepare_windowed first.
        if not prep_status.exists() or not window_path.exists():
            n_skip_noprep += 1
            sdir.mkdir(parents=True, exist_ok=True)
            sim_status.write_text(json.dumps({
                "scenario": sid, "ok": False,
                "skip_reason": "no_prep_outputs"}))
            print(f"[06sim] scenario {sid}: prep outputs missing; "
                  f"run prepare_windowed.py first", flush=True)
            continue
        try:
            prep = json.loads(prep_status.read_text())
            if not prep.get("ok"):
                n_skip_noprep += 1
                sim_status.write_text(json.dumps({
                    "scenario": sid, "ok": False,
                    "skip_reason": "prep_not_ok",
                    "prep_skip": prep.get("skip_reason"),
                    "prep_error": prep.get("error")}))
                continue
            window = json.loads(window_path.read_text())
            sim_start = window["sim_start"]
            sim_end = window["sim_end"]
            ftype = window["flow_type"]

            # Re-register the flow type for this process (the registry
            # is in-memory only; the on-disk inputs persist).
            register_flow_type(ftype, sin)

            import pywrdrb
            from src.pywrdrb.bridge import (
                SAVE_RESULTS_SETS, _get_parameter_subset_to_export,
            )
            simdir = sdir / "simulations"
            simdir.mkdir(parents=True, exist_ok=True)
            mb = pywrdrb.ModelBuilder(
                inflow_type=ftype, start_date=sim_start, end_date=sim_end,
                options={"inflow_ensemble_indices": ["0"],
                         "nyc_nj_demand_source": args.demand_source,
                         "flow_prediction_mode": "perfect_foresight"})
            mb.make_model()
            mfile = str(simdir / "model.json")
            mb.write_model(mfile)
            model = pywrdrb.Model.load(mfile)
            # Mirror run_pywrdrb_batch's exact output API: OutputRecorder
            # over the SAVE_RESULTS_SETS parameter subset BEFORE run().
            # Name batch_0.hdf5 so the aggregator's compute_metric_bank +
            # extract_drought_levels (glob batch_*.hdf5) reuse unchanged.
            opath = simdir / "batch_0.hdf5"
            all_param_names = [p.name for p in model.parameters if p.name]
            subset_names = _get_parameter_subset_to_export(
                all_param_names, SAVE_RESULTS_SETS)
            export_params = [p for p in model.parameters
                             if p.name in subset_names]
            recorder = pywrdrb.OutputRecorder(  # noqa: F841
                model=model, output_filename=str(opath),
                parameters=export_params)
            model.run()
            del recorder, model, mb

            sim_status.write_text(json.dumps({
                "scenario": sid, "ok": True,
                "sim_start": sim_start, "sim_end": sim_end,
                "demand_source": args.demand_source,
                "elapsed_s": round(time.time() - t0, 1)}))
            n_ok += 1
            print(f"[06sim] scenario {sid} OK ({sim_start}..{sim_end}) "
                  f"in {time.time()-t0:.1f}s", flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            n_fail += 1
            sdir.mkdir(parents=True, exist_ok=True)
            sim_status.write_text(json.dumps({
                "scenario": sid, "ok": False, "error": repr(e)}))
            print(f"[06sim] scenario {sid} FAILED: {e}\n"
                  f"{traceback.format_exc()}", flush=True)

    print(f"[06sim] array={args.array_id} rank={rank}: done "
          f"ok={n_ok} skip_done={n_skip_done} skip_noprep={n_skip_noprep} "
          f"fail={n_fail}", flush=True)


if __name__ == "__main__":
    main()
