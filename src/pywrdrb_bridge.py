"""Bridge between MOEA-FIND Pareto outputs and Pywr-DRB simulation inputs.

Stages 1-3 of the policy re-evaluation pipeline:

Stage 1 (generate):
    Replay Pareto DVs through multi-site Kirsch → Nowak daily → KDE downstream
    → marginal catchment inflows → node-first HDF5.

Stage 2 (prep):
    Run pywrdrb PredictedInflowEnsemblePreprocessor.

Stage 3 (simulate):
    Build and run Pywr-DRB model in batches, extract results.
"""

from __future__ import annotations

import gc
import json
import math
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from synhydro.core.ensemble import Ensemble
from synhydro.methods.disaggregation.temporal.nowak import NowakDisaggregator

from pywrdrb.pywr_drb_node_data import (
    immediate_downstream_nodes_dict,
    downstream_node_lags,
)
from pywrdrb.pre.flows import _subtract_upstream_catchment_inflows

from src.kirsch_wrapper import KirschBorgWrapper


# ===================================================================
# Stage 1: Generate multi-site daily traces from Pareto DVs
# ===================================================================


def replay_pareto_to_multisite_monthly(
    pareto_dvs: np.ndarray,
    wrapper: KirschBorgWrapper,
    kirsch_sites: List[str],
    start_date: str = "2030-01-01",
) -> Dict[int, pd.DataFrame]:
    """Replay each Pareto DV vector through multi-site Kirsch wrapper.

    Args:
        pareto_dvs: Array of shape ``(n_pareto, n_dvs)`` from Script 04.
        wrapper: :class:`KirschBorgWrapper` fitted on multi-site data.
        kirsch_sites: Ordered site names matching the generator's columns.
        start_date: Start date for the monthly index (calendar-year).

    Returns:
        ``{pareto_idx: DataFrame}`` with monthly datetime index,
        *kirsch_sites* columns, in MGD.
    """
    n_pareto = len(pareto_dvs)
    n_months = wrapper.n_years_out * 12
    dates = pd.date_range(start=start_date, periods=n_months, freq="MS")

    import time as _time
    monthly_traces: Dict[int, pd.DataFrame] = {}
    t0 = _time.time()
    log_every = max(1, n_pareto // 20)
    for i in range(n_pareto):
        arr = wrapper.generate(pareto_dvs[i])
        # Multi-site output: (n_years*12, n_sites), calendar-year order
        if arr.ndim == 1:
            arr = arr.reshape(n_months, -1)
        monthly_traces[i] = pd.DataFrame(arr, index=dates, columns=kirsch_sites)
        if (i + 1) % log_every == 0 or (i + 1) == n_pareto:
            elapsed = _time.time() - t0
            rate = (i + 1) / max(elapsed, 1e-6)
            eta = (n_pareto - (i + 1)) / max(rate, 1e-6)
            print(f"[bridge]   replay {i + 1}/{n_pareto} "
                  f"({rate:.1f} traces/s, ETA {eta:.0f}s)", flush=True)

    print(f"[bridge] Replayed {n_pareto} Pareto DVs → monthly multi-site "
          f"({n_months} months, {len(kirsch_sites)} sites) "
          f"in {_time.time() - t0:.1f}s", flush=True)
    return monthly_traces


def _disaggregate_chunk(
    chunk_traces: Dict[int, pd.DataFrame],
    nowak_disagg: NowakDisaggregator,
    seed: int,
) -> Dict[int, pd.DataFrame]:
    """Disaggregate one chunk of realizations inside a worker process.

    SynHydro's ``NowakDisaggregator.disaggregate`` hardcodes
    ``daily_realization_dict[0]`` when building its output metadata, so
    chunks whose first realization id is not 0 crash with
    ``KeyError: 0``. We work around this by re-keying the chunk to
    local 0..N-1 ids before calling SynHydro, then remapping the output
    back to the original global ids.
    """
    from synhydro.core.ensemble import EnsembleMetadata

    np.random.seed(seed)
    # Build a local-0-indexed view of the chunk, remembering the map
    # back to the caller's global ids.
    global_ids = list(chunk_traces.keys())
    local_chunk = {i: df for i, df in enumerate(chunk_traces.values())}
    sample_df = next(iter(local_chunk.values()))
    metadata = EnsembleMetadata(
        n_realizations=len(local_chunk),
        n_sites=len(sample_df.columns),
        time_resolution="MS",
    )
    monthly_ensemble = Ensemble(local_chunk, metadata=metadata)
    daily_ensemble = nowak_disagg.disaggregate(monthly_ensemble, seed=seed)
    # Remap local ids back to global ids
    out: Dict[int, pd.DataFrame] = {}
    for local_id, df in daily_ensemble.data_by_realization.items():
        out[int(global_ids[int(local_id)])] = df
    return out


def disaggregate_monthly_to_daily(
    monthly_traces: Dict[int, pd.DataFrame],
    nowak_disagg: NowakDisaggregator,
    seed: int = 42,
    n_jobs: int = -1,
    chunk_size: int = 50,
) -> Dict[int, pd.DataFrame]:
    """Nowak temporal disaggregation: monthly → daily for all Pareto traces.

    SynHydro's :meth:`NowakDisaggregator.disaggregate` iterates
    realizations serially inside a single process, which is the Stage 1
    bottleneck for large ensembles. We wrap it here with a joblib
    chunked parallel map so that the 16 allocated SLURM cores actually
    participate. Each chunk is an independent Ensemble disaggregated in
    its own worker process; results are merged on the caller side.

    Args:
        monthly_traces: ``{pareto_idx: monthly_DataFrame}`` from
            :func:`replay_pareto_to_multisite_monthly`.
        nowak_disagg: Fitted :class:`NowakDisaggregator`.
        seed: Base random seed. Per-chunk seeds are derived from this
            to keep reproducibility under parallel execution.
        n_jobs: Joblib worker count. ``-1`` uses all allocated CPUs
            (``SLURM_CPUS_PER_TASK``/``SLURM_NTASKS`` when set, otherwise
            ``os.cpu_count()``).
        chunk_size: Realizations per chunk. Smaller = finer progress
            granularity but more process-boundary overhead; 50 is a
            good balance at 7300 days × 22 sites.

    Returns:
        ``{pareto_idx: daily_DataFrame}`` with daily datetime index,
        same site columns, MGD units.
    """
    import os
    import time as _time

    from joblib import Parallel, delayed

    n_pareto = len(monthly_traces)
    if n_pareto == 0:
        return {}

    if n_jobs == -1:
        n_jobs = int(
            os.environ.get("SLURM_CPUS_PER_TASK")
            or os.environ.get("SLURM_NTASKS")
            or os.cpu_count() or 1
        )
    n_jobs = max(1, min(n_jobs, n_pareto))

    # Split into chunks keyed by original realization id
    items = sorted(monthly_traces.items(), key=lambda kv: int(kv[0]))
    chunks: List[Dict[int, pd.DataFrame]] = []
    for i in range(0, len(items), chunk_size):
        chunks.append(dict(items[i : i + chunk_size]))
    n_chunks = len(chunks)

    print(f"[bridge] Nowak disaggregation: {n_pareto} traces in "
          f"{n_chunks} chunks of up to {chunk_size}, n_jobs={n_jobs}",
          flush=True)

    t0 = _time.time()
    # ``max_nbytes=None`` disables joblib's automatic read-only memory
    # mapping of large arrays. SynHydro's NowakDisaggregator does
    # in-place arithmetic on its internal arrays (e.g.
    # ``site_props[mask] *= ...``), which raises
    # ``ValueError: assignment destination is read-only`` when the
    # arrays arrive memory-mapped in a worker. The extra per-worker
    # memory cost is acceptable at this ensemble size.
    results = Parallel(
        n_jobs=n_jobs, backend="loky", verbose=0, max_nbytes=None,
    )(
        delayed(_disaggregate_chunk)(chunk, nowak_disagg, seed + k)
        for k, chunk in enumerate(chunks)
    )
    elapsed = _time.time() - t0

    # Merge chunk outputs
    daily_traces: Dict[int, pd.DataFrame] = {}
    for r in results:
        daily_traces.update(r)

    sample_df = next(iter(daily_traces.values()))
    print(f"[bridge] Nowak disaggregation complete: "
          f"{len(daily_traces)} traces, {len(sample_df)} days, "
          f"{len(sample_df.columns)} sites in {elapsed:.1f}s "
          f"({len(daily_traces) / max(elapsed, 1e-6):.1f} traces/s)",
          flush=True)
    return daily_traces


def fit_kde_models(
    Q_inflow_hist: pd.DataFrame,
    kde_pairs: List[Tuple[str, str]],
) -> Dict[str, stats.gaussian_kde]:
    """Fit KDE models on historical upstream/downstream flow ratios.

    For each ``(upstream, downstream)`` pair, fits a Gaussian KDE on
    the ratio ``downstream_inflow / upstream_flow`` from the historical
    catchment inflow data. Follows the pattern from
    ``StochasticExploratoryExperiment/methods/generate.py``.

    Args:
        Q_inflow_hist: Historical catchment inflow DataFrame (daily, MGD).
        kde_pairs: List of ``(upstream, downstream)`` node name pairs.

    Returns:
        ``{f"{upstream}_to_{downstream}": fitted_kde}`` dict.
    """
    kdes: Dict[str, stats.gaussian_kde] = {}
    for upstream, downstream in kde_pairs:
        if upstream not in Q_inflow_hist.columns:
            print(f"[bridge] WARNING: upstream {upstream} not in inflow data, "
                  f"skipping KDE for {upstream}→{downstream}")
            continue
        if downstream not in Q_inflow_hist.columns:
            print(f"[bridge] WARNING: downstream {downstream} not in inflow data, "
                  f"skipping KDE for {upstream}→{downstream}")
            continue

        xs = Q_inflow_hist[upstream]
        ys = Q_inflow_hist[downstream]
        frac = ys / xs
        frac = frac[~np.isnan(frac)]
        frac = frac[np.isfinite(frac)]

        if len(frac) < 10:
            print(f"[bridge] WARNING: only {len(frac)} valid ratios for "
                  f"{upstream}→{downstream}, skipping")
            continue

        kde = stats.gaussian_kde(frac.values)
        kde_name = f"{upstream}_to_{downstream}"
        kdes[kde_name] = kde

    print(f"[bridge] Fitted {len(kdes)} KDE models for downstream regression")
    return kdes


def generate_kde_downstream_nodes(
    daily_traces: Dict[int, pd.DataFrame],
    kdes: Dict[str, stats.gaussian_kde],
    kde_pairs: List[Tuple[str, str]],
    seed: int = 42,
) -> Dict[int, pd.DataFrame]:
    """Generate flows at non-major nodes via KDE regression.

    For each ``(upstream, downstream)`` pair:
        1. Sample ratios from the fitted KDE
        2. ``downstream_inflow = upstream_flow * ratio``
        3. ``downstream_gage_flow = downstream_inflow + lagged_upstream_flow``

    Modifies *daily_traces* in place (adds downstream columns) and returns
    the same dict.

    Args:
        daily_traces: ``{pareto_idx: daily_DataFrame}`` with Kirsch-generated
            sites as columns.
        kdes: Fitted KDE models from :func:`fit_kde_models`.
        kde_pairs: ``(upstream, downstream)`` pairs to generate.
        seed: Base seed for reproducibility (varied per pair and solution).

    Returns:
        ``{pareto_idx: daily_DataFrame}`` with both Kirsch-generated and
        KDE-regressed columns.
    """
    import time as _time
    n_pareto = len(daily_traces)
    sample_df = next(iter(daily_traces.values()))
    n_days = len(sample_df)
    t0 = _time.time()

    for pair_idx, (upstream, downstream) in enumerate(kde_pairs):
        kde_name = f"{upstream}_to_{downstream}"
        if kde_name not in kdes:
            continue

        kde = kdes[kde_name]
        lag = downstream_node_lags.get(downstream, 0)

        # Draw all samples at once: (n_days, n_pareto)
        kde_seed = seed + hash(kde_name) % 10000
        n_total = n_days * n_pareto
        samples = kde.resample(n_total, seed=kde_seed).flatten()
        samples = np.clip(samples, 0, 1)
        samples = samples.reshape(n_days, n_pareto)

        for col_idx, (pareto_idx, df) in enumerate(daily_traces.items()):
            if upstream not in df.columns:
                continue
            upstream_flow = df[upstream].values
            ratios = samples[:, col_idx]
            downstream_inflow = upstream_flow * ratios

            # Add lagged upstream contribution for gage flow
            if lag > 0:
                downstream_gage = downstream_inflow.copy()
                downstream_gage[lag:] += upstream_flow[:-lag]
                downstream_gage[:lag] += upstream_flow[:lag]
            else:
                downstream_gage = downstream_inflow + upstream_flow

            df[downstream] = downstream_gage
        print(f"[bridge]   KDE pair {pair_idx + 1}/{len(kde_pairs)} "
              f"({upstream}→{downstream}) done "
              f"at {_time.time() - t0:.1f}s", flush=True)

    # Add delTrenton column (zeroed by convention)
    for df in daily_traces.values():
        df["delTrenton"] = 0.0

    print(f"[bridge] KDE regression added {len(kde_pairs)} downstream nodes "
          f"+ delTrenton in {_time.time() - t0:.1f}s",
          flush=True)
    return daily_traces


def compute_marginal_catchment_inflows(
    gage_flow_traces: Dict[int, pd.DataFrame],
) -> Dict[int, pd.DataFrame]:
    """Convert gage flows to marginal catchment inflows.

    Uses ``pywrdrb.pre.flows._subtract_upstream_catchment_inflows()``
    which subtracts upstream contributions (with lag) from each node's
    cumulative gage flow, yielding the direct catchment inflow at each node.

    Args:
        gage_flow_traces: ``{pareto_idx: gage_flow_DataFrame}``.

    Returns:
        ``{pareto_idx: catchment_inflow_DataFrame}`` with marginal inflows.
    """
    import time as _time
    t0 = _time.time()
    inflow_traces: Dict[int, pd.DataFrame] = {}
    n = len(gage_flow_traces)
    log_every = max(1, n // 10)
    for i, (idx, gage_df) in enumerate(gage_flow_traces.items()):
        inflow_traces[idx] = _subtract_upstream_catchment_inflows(gage_df)
        if (i + 1) % log_every == 0 or (i + 1) == n:
            print(f"[bridge]   catchment-inflow {i + 1}/{n} "
                  f"at {_time.time() - t0:.1f}s", flush=True)

    print(f"[bridge] Computed marginal catchment inflows for "
          f"{len(inflow_traces)} traces in {_time.time() - t0:.1f}s",
          flush=True)
    return inflow_traces


def write_flowensemble_hdf5(
    traces: Dict[int, pd.DataFrame],
    output_path: Path,
) -> None:
    """Write traces as node-first HDF5 matching Pywr-DRB FlowEnsemble format.

    Uses SynHydro's :class:`Ensemble` with ``stored_by_node=True`` which
    produces the exact HDF5 structure that
    ``pywrdrb.parameters.ensemble.FlowEnsemble`` expects:

    .. code-block::

        /cannonsville/
            attrs["column_labels"] = ["0", "1", ...]
            dataset "date": ISO date strings
            dataset "0": float64 array
            dataset "1": float64 array
        /pepacton/
            ...

    Args:
        traces: ``{pareto_idx: DataFrame}`` — each DataFrame has site name
            columns and a daily datetime index.
        output_path: Path to output HDF5 file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ensemble = Ensemble(traces)
    ensemble.to_hdf5(str(output_path), stored_by_node=True)

    print(f"[bridge] Wrote {output_path.name} "
          f"({len(traces)} realizations, "
          f"{len(next(iter(traces.values())).columns)} nodes)")


# ===================================================================
# Stage 2: Prepare Pywr-DRB preprocessor inputs
# ===================================================================


def register_flow_type(flow_type: str, ensemble_dir: Path) -> None:
    """Register a custom flow type directory with pywrdrb's path navigator.

    After registration, ``pn.sc.get(f"flows/{flow_type}")`` resolves to
    *ensemble_dir*, allowing ``ModelBuilder(inflow_type=flow_type)`` to
    find the HDF5 files.

    Args:
        flow_type: Name for this flow dataset (e.g., ``"moea_find_pareto"``).
        ensemble_dir: Directory containing ``catchment_inflow_mgd.hdf5`` and
            related files.
    """
    import pywrdrb

    pn_config = pywrdrb.get_pn_config()
    pn_config[f"flows/{flow_type}"] = str(Path(ensemble_dir).resolve())
    pywrdrb.load_pn_config(pn_config)

    print(f"[bridge] Registered flow type '{flow_type}' → {ensemble_dir}")


def _get_mpi_context():
    """Return ``(comm, rank, size)`` if MPI is available, else serial fallback."""
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        return comm, comm.Get_rank(), comm.Get_size()
    except ImportError:
        return None, 0, 1


def prep_predicted_inflows(
    flow_type: str,
    ensemble_dir: Path,
    realization_ids: List[str],
    use_mpi: bool = False,
    demand_source: str = "constant_max",
) -> None:
    """Run Pywr-DRB ensemble preprocessors needed for the upcoming sim.

    Which files ModelBuilder actually reads depends on the demand source:

    +----------------------------------------+---------------+----------+
    | Artifact                               | constant_max  | custom / |
    |                                        |               | historic |
    +----------------------------------------+---------------+----------+
    | predicted_inflows_mgd.hdf5             | required      | required |
    | diversion_nj_extrapolated_mgd.hdf5     | unused        | required |
    | diversion_nyc_extrapolated_mgd.hdf5    | unused        | required |
    | predicted_diversions_mgd.hdf5          | unused[1]     | required |
    +----------------------------------------+---------------+----------+

    [1] Requires the Pywr-DRB ModelBuilder patch that substitutes a
    ``constant`` parameter for ``predicted_demand_nj_lag{1..4}`` when
    ``nyc_nj_demand_source == "constant_max"``. With that patch in
    place, the three diversion HDF5s are not needed under constant_max
    and we skip them here. Without the patch, ModelBuilder would still
    demand them even though they encode a forecast that is
    operationally unused; force-regenerate by passing a non-constant
    ``demand_source``.

    Args:
        flow_type: Registered flow type name.
        ensemble_dir: Directory containing ``catchment_inflow_mgd.hdf5``
            and ``gage_flow_mgd.hdf5`` (both written by Stage 1).
        realization_ids: List of realization ID strings (e.g.,
            ``["0","1",...]``).
        use_mpi: If True, distribute preprocessing across MPI ranks.
        demand_source: Selects which preprocessors run. ``constant_max``
            runs only Step 1; any other value runs all four steps to
            match ModelBuilder's `historical` / `custom` demand paths.
    """
    from pywrdrb.pre import (
        PredictedInflowEnsemblePreprocessor,
        ExtrapolatedDiversionEnsemblePreprocessor,
        PredictedDiversionEnsemblePreprocessor,
    )

    comm = None
    if use_mpi:
        comm, rank, size = _get_mpi_context()
        use_mpi = size > 1
    else:
        rank = 0

    catchment_inflow_file = str(Path(ensemble_dir) / "catchment_inflow_mgd.hdf5")
    gage_flow_file = str(Path(ensemble_dir) / "gage_flow_mgd.hdf5")
    run_diversion_steps = demand_source != "constant_max"
    n_steps = 4 if run_diversion_steps else 1

    if rank == 0:
        print(f"[bridge] Step 1/{n_steps}: PredictedInflowEnsemblePreprocessor "
              f"for {len(realization_ids)} realizations "
              f"({'MPI' if use_mpi else 'serial'})")
    inflow_pre = PredictedInflowEnsemblePreprocessor(
        flow_type=flow_type,
        ensemble_hdf5_file=catchment_inflow_file,
        realization_ids=realization_ids,
        start_date=None, end_date=None,
        modes=("perfect_foresight",),
        use_log=True, remove_zeros=True, use_const=False,
        use_mpi=use_mpi, comm=comm,
    )
    inflow_pre.load(); inflow_pre.process(); inflow_pre.save()
    del inflow_pre

    if not run_diversion_steps:
        if rank == 0:
            print(f"[bridge] demand_source={demand_source!r}: skipping "
                  f"diversion preprocessors (ModelBuilder uses constant "
                  f"scalars for demand_nj/demand_nyc and predicted_demand_nj_lag*).")
            print(f"[bridge] Wrote predicted_inflows_mgd.hdf5 → {ensemble_dir}")
        return

    if rank == 0:
        print(f"[bridge] Step 2/{n_steps}: ExtrapolatedDiversionEnsemblePreprocessor(nj)")
    nj_extrap = ExtrapolatedDiversionEnsemblePreprocessor(
        loc="nj",
        flow_type=flow_type,
        ensemble_hdf5_file=gage_flow_file,
        realization_ids=realization_ids,
        use_mpi=use_mpi, comm=comm,
    )
    nj_extrap.load(); nj_extrap.process(); nj_extrap.save()
    del nj_extrap

    if rank == 0:
        print(f"[bridge] Step 3/{n_steps}: ExtrapolatedDiversionEnsemblePreprocessor(nyc)")
    nyc_extrap = ExtrapolatedDiversionEnsemblePreprocessor(
        loc="nyc",
        flow_type=flow_type,
        ensemble_hdf5_file=gage_flow_file,
        realization_ids=realization_ids,
        use_mpi=use_mpi, comm=comm,
    )
    nyc_extrap.load(); nyc_extrap.process(); nyc_extrap.save()
    del nyc_extrap

    if rank == 0:
        print(f"[bridge] Step 4/{n_steps}: PredictedDiversionEnsemblePreprocessor")
    import pywrdrb
    nj_div_hdf5 = str(
        pywrdrb.get_pn_object().sc.get(f"flows/{flow_type}")
        / "diversion_nj_extrapolated_mgd.hdf5"
    )
    div_pre = PredictedDiversionEnsemblePreprocessor(
        flow_type=flow_type,
        ensemble_hdf5_file=nj_div_hdf5,
        realization_ids=realization_ids,
        start_date=None, end_date=None,
        modes=("perfect_foresight",),
        use_log=True, remove_zeros=True, use_const=False,
        use_mpi=use_mpi, comm=comm,
    )
    div_pre.load(); div_pre.process(); div_pre.save()
    del div_pre

    if rank == 0:
        print(f"[bridge] Wrote predicted_inflows_mgd.hdf5, "
              f"diversion_nj_extrapolated_mgd.hdf5, "
              f"diversion_nyc_extrapolated_mgd.hdf5, "
              f"predicted_diversions_mgd.hdf5 → {ensemble_dir}")


# ===================================================================
# Stage 3: Run Pywr-DRB simulations
# ===================================================================

# Results sets to record. Keep aligned with
# ``../StochasticExploratoryExperiment/methods/config.py::SAVE_RESULTS_SETS``
# so Stage 4 metric-bank helpers can compute FFMP exposure, Hashimoto
# reliability/vulnerability, flow-target reliability, and any storage- or
# delivery-based satisficing metric from a single pass over the HDF5
# output. Extending this list only adds columns; it does not change the
# Pywr-DRB model or the simulation cost.
SAVE_RESULTS_SETS = [
    "major_flow",
    "inflow",
    "res_storage",
    "lower_basin_mrf_contributions",
    "mrf_target",
    "ibt_diversions",
    "ibt_demands",
    "nyc_release_components",
    "res_level",
]


def _get_parameter_subset_to_export(all_parameter_names, results_sets):
    """Filter parameter names to those belonging to *results_sets*.

    Mirrors ``StochasticExploratoryExperiment/methods/utils.py::
    get_parameter_subset_to_export``.
    """
    import pywrdrb
    output_loader = pywrdrb.load.Output(output_filenames=[])
    keep_keys = []
    for rs in results_sets:
        if rs == "all":
            continue
        keys_sub, _ = output_loader.get_keys_and_column_names_for_results_set(
            all_parameter_names, rs,
        )
        keep_keys.extend(keys_sub)
    return keep_keys


def run_pywrdrb_batch(
    flow_type: str,
    realization_ids: List[str],
    start_date: str,
    end_date: str,
    output_dir: Path,
    model_dir: Path,
    batch_size: int = 10,
    demand_source: str = "constant_max",
    use_mpi: bool = False,
    combined_output_name: str = "pywrdrb_output.hdf5",
) -> Path:
    """Run Pywr-DRB simulations in batches, optionally parallelized via MPI.

    When ``use_mpi=True``, realizations are distributed across MPI ranks;
    each rank runs its subset in serial batches of *batch_size*.  After all
    ranks finish, rank 0 combines the batch HDF5 files.

    When ``use_mpi=False`` (default), everything runs on a single process.

    Args:
        flow_type: Registered flow type name.
        realization_ids: All realization ID strings.
        start_date: Pywr-DRB simulation start (e.g., ``"2030-01-01"``).
        end_date: Pywr-DRB simulation end (e.g., ``"2049-12-31"``).
        output_dir: Directory for simulation output HDF5 files.
        model_dir: Directory for temporary model JSON files.
        batch_size: Number of realizations per Pywr-DRB run.
        demand_source: NYC/NJ demand source option (default ``"constant_max"``).
        use_mpi: Distribute across MPI ranks.
        combined_output_name: Filename for the combined output HDF5.

    Returns:
        Path to the combined output HDF5 file.
    """
    import glob
    import re
    import pywrdrb
    from pywrdrb.utils.hdf5 import combine_batched_hdf5_outputs

    output_dir = Path(output_dir)
    model_dir = Path(model_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    model_dir.mkdir(parents=True, exist_ok=True)

    combined_output = output_dir / combined_output_name

    # --- MPI context ---
    comm = None
    if use_mpi:
        comm, rank, size = _get_mpi_context()
        use_mpi = size > 1
    else:
        rank, size = 0, 1

    # --- Distribute realizations across ranks ---
    rank_realization_ids = np.array_split(realization_ids, size)[rank]
    rank_realization_ids = list(rank_realization_ids)

    if rank == 0:
        print(f"[bridge] Running Pywr-DRB: {len(realization_ids)} realizations, "
              f"{size} rank(s), batch_size={batch_size}")

    # --- Clean old batch files (rank 0 only) ---
    if rank == 0:
        for old in glob.glob(str(output_dir / f"rank*_batch*.hdf5")):
            os.remove(old)
        for old in glob.glob(str(model_dir / f"model_rank*_batch*.json")):
            os.remove(old)

    if comm is not None:
        comm.Barrier()

    # --- Run batches for this rank's realizations ---
    n_rank_real = len(rank_realization_ids)
    n_batches = math.ceil(n_rank_real / batch_size) if n_rank_real > 0 else 0

    for batch_idx in range(n_batches):
        start = batch_idx * batch_size
        end_idx = min(start + batch_size, n_rank_real)
        batch_ids = rank_realization_ids[start:end_idx]

        print(f"[bridge] rank={rank} batch {batch_idx + 1}/{n_batches}: "
              f"realizations {batch_ids[0]}–{batch_ids[-1]}")

        model_options = {
            "inflow_ensemble_indices": batch_ids,
            "nyc_nj_demand_source": demand_source,
            "flow_prediction_mode": "perfect_foresight",
        }

        mb = pywrdrb.ModelBuilder(
            inflow_type=flow_type,
            start_date=start_date,
            end_date=end_date,
            options=model_options,
        )
        mb.make_model()

        model_fname = str(model_dir / f"model_rank{rank}_batch{batch_idx}.json")
        mb.write_model(model_fname)
        model = pywrdrb.Model.load(model_fname)

        # Filter to results sets we care about
        all_param_names = [p.name for p in model.parameters if p.name]
        subset_names = _get_parameter_subset_to_export(
            all_param_names, SAVE_RESULTS_SETS,
        )
        export_params = [p for p in model.parameters
                         if p.name in subset_names]

        batch_output = str(
            output_dir / f"rank{rank}_batch{batch_idx}.hdf5"
        )
        recorder = pywrdrb.OutputRecorder(
            model=model,
            output_filename=batch_output,
            parameters=export_params,
        )

        model.run()

        # Drop every object that pins the Model/Recorder graph before the
        # next batch's ModelBuilder allocates. Without this, rank-level
        # memory grew batch-over-batch and OOM-killed one rank on 212267
        # partway through batch 2/2.
        del recorder, model, mb, export_params, all_param_names
        gc.collect()

        print(f"[bridge]   rank={rank} batch {batch_idx + 1} complete")

    # --- Wait for all ranks ---
    if comm is not None:
        comm.Barrier()

    # --- Combine on rank 0 ---
    if rank == 0:
        batch_pattern = str(output_dir / "rank*_batch*.hdf5")
        all_batch_files = sorted(
            glob.glob(batch_pattern),
            key=lambda f: tuple(
                int(x) for x in re.findall(r'\d+', Path(f).stem)
            ),
        )
        if not all_batch_files:
            print("[bridge] WARNING: no batch output files found!")
            return combined_output

        print(f"[bridge] Combining {len(all_batch_files)} batch files...")

        if combined_output.exists():
            combined_output.unlink()
        combine_batched_hdf5_outputs(all_batch_files, str(combined_output))

        # Cleanup batch + model files
        for f in all_batch_files:
            os.remove(f)
        for f in glob.glob(str(model_dir / "model_rank*_batch*.json")):
            os.remove(f)

        print(f"[bridge] Combined output: {combined_output}")

    return combined_output
