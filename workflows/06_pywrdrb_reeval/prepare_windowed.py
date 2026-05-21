"""Stage 06 (windowed) -- per-scenario Pywr-DRB INPUT PREPARATION.

Split out of ``policy_reeval_windowed.py`` (2026-05-20) so the heavy
Kirsch-replay / Nowak-disaggregation / KDE chain is cached on disk and
the simulation step can be re-run without redoing it. Pair with:

    1) prepare_windowed.py    (this file)  -- DV -> daily multisite -> HDF5
    2) policy_reeval_windowed.py           -- HDF5 + window -> pywrdrb sim
    3) aggregate_windowed.py               -- per-scenario -> metric bank

Per scenario this script does, in one pass:

  replay DV -> multisite monthly -> SSI-3 (ref site) -> first-event
  window (sim_start, sim_end) -> daily disaggregation (Nowak) -> KDE
  downstream nodes -> slice every site to the window -> compute
  marginal catchment inflows -> write gage_flow_mgd.hdf5 +
  catchment_inflow_mgd.hdf5 -> register_flow_type + prep_predicted_inflows

Outputs (per scenario ``s`` under the upstream MOEA slug ``src_slug``):
  outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/
    scenarios/<s>/
      pywrdrb_inputs/{gage_flow_mgd.hdf5, catchment_inflow_mgd.hdf5,
                      predicted_inflows/...}
      window.json     (scenario, src_slug, sim_start, sim_end, window meta)
      status_prep.json (ok / skipped / failed)

Design: EMBARRASSINGLY PARALLEL, no MPI collectives. Work is split by
(SLURM array task) x (MPI rank, 1 core each); each rank owns a disjoint
scenario subset. Mirrors ``policy_reeval_windowed.py`` topology so the
two stages can re-use the same SLURM allocation pattern.
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
    ap.add_argument("--demand-source", default="constant_max",
                    help="Demand label baked into the predicted_inflows "
                         "preprocessor (the simulation step can also "
                         "override --demand-source at run time).")
    ap.add_argument("--disagg-seed", type=int, default=42)
    ap.add_argument("--array-id", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--array-count", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1)))
    ap.add_argument("--use-mpi", action="store_true")
    ap.add_argument("--overwrite", action="store_true",
                    help="Re-prep scenarios whose status_prep.json reports "
                         "ok; default behavior skips them.")
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
    print(f"[06prep] array={args.array_id}/{args.array_count} "
          f"rank={rank}/{size} -> {len(my_ids)} scenarios of {n_total} "
          f"(src={src_slug})", flush=True)
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
    # stage-04 optimised, else it is inconsistent with the objectives.
    # Reproduce stage-04's exact path:
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
    n_ok, n_skip_done, n_skip_empty, n_fail = 0, 0, 0, 0
    for sid in my_ids:
        t0 = time.time()
        sdir = scen_root / str(sid)
        sin = sdir / "pywrdrb_inputs"
        prep_status = sdir / "status_prep.json"

        if not args.overwrite and prep_status.exists():
            try:
                if json.loads(prep_status.read_text()).get("ok"):
                    n_skip_done += 1
                    continue
            except Exception:  # noqa: BLE001
                pass  # malformed status -> re-prep

        try:
            dv = pareto_dvs[sid:sid + 1]
            # Window from the MOEA single-site SSI-3 (stage-04-consistent).
            syn_moea = np.asarray(
                moea_gen.generate(pareto_dvs[sid]), dtype=float).flatten()
            ssi = ssi_calc.transform(
                flows_to_series(syn_moea, start_date="2100-01-01"))
            win = first_event_window(ssi)
            if win is None:
                n_skip_empty += 1
                sdir.mkdir(parents=True, exist_ok=True)
                prep_status.write_text(json.dumps({
                    "scenario": sid, "ok": False,
                    "skip_reason": "no_critical_event",
                    "elapsed_s": round(time.time() - t0, 1)}))
                print(f"[06prep] scenario {sid}: no critical event; skipped",
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
            sin.mkdir(parents=True, exist_ok=True)
            write_flowensemble_hdf5(win_full, sin / "gage_flow_mgd.hdf5")
            write_flowensemble_hdf5(catch, sin / "catchment_inflow_mgd.hdf5")

            register_flow_type(ftype, sin)
            prep_predicted_inflows(ftype, sin, ["0"], use_mpi=False,
                                   demand_source=args.demand_source)

            (sdir / "window.json").write_text(json.dumps({
                "scenario": sid, "src_slug": src_slug,
                "sim_start": sim_start, "sim_end": sim_end,
                "flow_type": ftype,
                **win.to_dict()}, indent=2))
            prep_status.write_text(json.dumps({
                "scenario": sid, "ok": True,
                "flow_type": ftype,
                "sim_start": sim_start, "sim_end": sim_end,
                "elapsed_s": round(time.time() - t0, 1)}))
            n_ok += 1
            print(f"[06prep] scenario {sid} OK ({sim_start}..{sim_end}, "
                  f"{win.n_window_months}mo, "
                  f"spinup_short={win.short_spinup}) "
                  f"in {time.time()-t0:.1f}s", flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            n_fail += 1
            sdir.mkdir(parents=True, exist_ok=True)
            prep_status.write_text(json.dumps({
                "scenario": sid, "ok": False, "error": repr(e)}))
            print(f"[06prep] scenario {sid} FAILED: {e}\n"
                  f"{traceback.format_exc()}", flush=True)

    print(f"[06prep] array={args.array_id} rank={rank}: done "
          f"ok={n_ok} skip_done={n_skip_done} skip_empty={n_skip_empty} "
          f"fail={n_fail}", flush=True)


if __name__ == "__main__":
    main()
