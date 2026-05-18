"""Stage 06 (windowed) -- per-scenario first-drought-window Pywr-DRB re-eval.

Replaces the monolithic 4-stage policy_reeval for the windowed pipeline
(2026-05-15 redesign). Each Pareto scenario is trimmed to its first
critical SSI-3 drought window + buffer (src.metrics.drought_window) and run
through Pywr-DRB independently with its own start/end date and
batch_size=1.

Design: EMBARRASSINGLY PARALLEL, ZERO shared state / MPI collectives.
Work is split by (SLURM array task) x (MPI rank, 1 core each); each rank
owns a disjoint scenario subset and does, per scenario, the full chain:

  replay DV -> multisite monthly -> SSI-3 (ref site) -> first-event
  window -> daily disaggregation -> KDE downstream -> per-scenario
  pywrdrb inputs -> STARFIT presim + predicted inflows (serial, tiny) ->
  pywrdrb ModelBuilder run over the window -> per-scenario metric row.

Outputs (per scenario ``s`` under the COARSE src slug):
  outputs/06_pywrdrb_reeval/policy_reeval/<slug>/
    scenarios/<s>/{pywrdrb_inputs/, simulations/pywrdrb_output.hdf5,
                   metrics.json, window.json}
A separate aggregation pass (aggregate_windowed.py) concatenates the
per-scenario metrics + window manifests into the standard
results/metric_bank.parquet that stages 07/08/09 consume.

This script processes ONLY this worker's scenario subset and exits;
it never combines outputs (the aggregator does, after all array tasks).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.multisite_data import (  # noqa: E402
    load_pywrdrb_gage_flow,
    load_pywrdrb_catchment_inflow,
    get_kirsch_sites,
    get_kde_regression_sites,
    get_kde_pairs,
    fit_multisite_generators,
)
from src.pywrdrb.bridge import (  # noqa: E402
    fit_kde_models,
    replay_pareto_to_multisite_monthly,
    disaggregate_monthly_to_daily,
    generate_kde_downstream_nodes,
    compute_marginal_catchment_inflows,
    write_flowensemble_hdf5,
    prep_predicted_inflows,
    register_flow_type,
    _get_mpi_context,
)
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.metrics.objectives import make_ssi_calculator, flows_to_series  # noqa: E402
from src.metrics.drought_window import first_event_window, window_dates  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "06_pywrdrb_reeval"
DRIVER = "policy_reeval"


def _worker_scenarios(n_total: int, array_id: int, array_count: int,
                      rank: int, size: int) -> list[int]:
    """Disjoint scenario subset for (array task, MPI rank).

    First split 0..n_total across array tasks, then across ranks within
    the task. Pure index arithmetic -- no communication.
    """
    task_chunk = np.array_split(np.arange(n_total), array_count)[array_id]
    return [int(i) for i in np.array_split(task_chunk, size)[rank]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pareto-results", type=Path, required=True)
    ap.add_argument("--site-label", default="cannonsville",
                    help="Reference site whose SSI-3 defines the window "
                         "(must match the stage-04 objective site).")
    ap.add_argument("--ssi", type=int, default=3)
    ap.add_argument("--baseline-dataset", default="nhmv10")
    ap.add_argument("--demand-source", default="constant_max")
    ap.add_argument("--disagg-seed", type=int, default=42)
    ap.add_argument("--array-id", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--array-count", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1)))
    ap.add_argument("--use-mpi", action="store_true")
    args = ap.parse_args()

    comm = None
    rank, size = 0, 1
    if args.use_mpi:
        comm, rank, size = _get_mpi_context()

    results = json.loads(Path(args.pareto_results).read_text())
    pareto_dvs = np.asarray(results["pareto_dvs"], dtype=float)
    n_years = int(results["n_years_out"])
    n_total = len(pareto_dvs)
    src_slug = args.pareto_results.parent.name

    # results.json["mode"] is the generator *label* (e.g. "Kirsch
    # (residual)"), not the bare wrapper mode token. Prefer the clean
    # value from the sibling stage-04 config.json; else sanitise.
    cfg_path = Path(args.pareto_results).parent / "config.json"
    if cfg_path.exists():
        mode = json.loads(cfg_path.read_text()).get("mode", "residual")
    else:
        raw_mode = str(results.get("mode", "residual")).lower()
        mode = "index" if "index" in raw_mode else "residual"
    out = stage_output_dir(STAGE, DRIVER, src_slug)
    flow_type_base = f"moea_find_{src_slug}"

    my_ids = _worker_scenarios(n_total, args.array_id, args.array_count,
                               rank, size)
    print(f"[06w] array={args.array_id}/{args.array_count} rank={rank}/{size} "
          f"-> {len(my_ids)} scenarios of {n_total} (src={src_slug})",
          flush=True)
    if not my_ids:
        return

    # ---- Shared fits (per worker; read-only, no cross-rank comm) ----
    Q_gage = load_pywrdrb_gage_flow(args.baseline_dataset)
    Q_inflow = load_pywrdrb_catchment_inflow(args.baseline_dataset)
    kirsch_sites = get_kirsch_sites(Q_gage)
    kde_sites = get_kde_regression_sites(Q_gage)
    kde_pairs = get_kde_pairs(kirsch_sites, kde_sites)
    kirsch_gen, nowak_disagg = fit_multisite_generators(Q_gage, kirsch_sites)
    kdes = fit_kde_models(Q_inflow, kde_pairs)
    wrapper = KirschBorgWrapper(kirsch_gen, mode=mode, n_years_out=n_years)

    # The window MUST be defined by the SAME single-site SSI-3 first event
    # stage-04 optimised, else it is inconsistent with the objectives (and
    # the bridge's multisite/MGD replay finds no event at all -> every
    # scenario wrongly skipped). Reproduce stage-04's exact path:
    # prepare_data -> build_kirsch_generator -> KirschBorgWrapper (single
    # MOEA site) -> SSI-3 prefit on historical monthly_1d.
    from src.experiment import prepare_data
    from src.hydrology.kirsch_utils import build_kirsch_generator
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    moea_gen = KirschBorgWrapper(build_kirsch_generator(monthly_2d),
                                 mode=mode, n_years_out=n_years)
    ssi_calc = make_ssi_calculator(timescale=args.ssi)
    ssi_calc.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))

    scen_root = out / "scenarios"
    n_ok, n_skip = 0, 0
    for sid in my_ids:
        t0 = time.time()
        sdir = scen_root / str(sid)
        try:
            dv = pareto_dvs[sid:sid + 1]
            # Window from the MOEA single-site SSI-3 (stage-04-consistent).
            syn_moea = np.asarray(
                moea_gen.generate(pareto_dvs[sid]), dtype=float).flatten()
            ssi = ssi_calc.transform(
                flows_to_series(syn_moea, start_date="2100-01-01"))
            win = first_event_window(ssi)
            if win is None:
                n_skip += 1
                print(f"[06w] scenario {sid}: no critical event; skipped",
                      flush=True)
                continue
            dts = window_dates(win, trace_start="2030-01-01")
            sim_start = dts["sim_start"].strftime("%Y-%m-%d")
            sim_end = dts["sim_end"].strftime("%Y-%m-%d")

            # Multisite replay of the SAME DV for the pywrdrb daily inputs
            # (separate from the single-site MOEA generator used for the
            # window above; both consume pareto_dvs[sid]).
            monthly = replay_pareto_to_multisite_monthly(
                dv, wrapper, kirsch_sites, start_date="2030-01-01")
            daily = disaggregate_monthly_to_daily(
                monthly, nowak_disagg, seed=args.disagg_seed, n_jobs=1)
            daily_full = generate_kde_downstream_nodes(
                daily, kdes, kde_pairs, seed=args.disagg_seed)
            # Slice every site's daily series to the window.
            win_full = {0: daily_full[0].loc[sim_start:sim_end]}
            catch = compute_marginal_catchment_inflows(win_full)

            ftype = f"{flow_type_base}_s{sid}"
            sin = sdir / "pywrdrb_inputs"
            sin.mkdir(parents=True, exist_ok=True)
            write_flowensemble_hdf5(win_full, sin / "gage_flow_mgd.hdf5")
            write_flowensemble_hdf5(catch, sin / "catchment_inflow_mgd.hdf5")

            from src.pywrdrb.bridge import register_flow_type
            register_flow_type(ftype, sin)
            prep_predicted_inflows(ftype, sin, ["0"], use_mpi=False,
                                   demand_source=args.demand_source)

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

            (sdir / "window.json").write_text(json.dumps({
                "scenario": sid, "src_slug": src_slug,
                "sim_start": sim_start, "sim_end": sim_end,
                **win.to_dict()}, indent=2))
            (sdir / "status.json").write_text(json.dumps({
                "scenario": sid, "ok": True,
                "elapsed_s": round(time.time() - t0, 1)}))
            n_ok += 1
            print(f"[06w] scenario {sid} OK ({sim_start}..{sim_end}, "
                  f"{win.n_window_months}mo, spinup_short={win.short_spinup}) "
                  f"in {time.time()-t0:.1f}s", flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            n_skip += 1
            sdir.mkdir(parents=True, exist_ok=True)
            (sdir / "status.json").write_text(json.dumps({
                "scenario": sid, "ok": False, "error": repr(e)}))
            print(f"[06w] scenario {sid} FAILED: {e}\n"
                  f"{traceback.format_exc()}", flush=True)

    print(f"[06w] array={args.array_id} rank={rank}: done "
          f"ok={n_ok} skip/fail={n_skip}", flush=True)


if __name__ == "__main__":
    main()
