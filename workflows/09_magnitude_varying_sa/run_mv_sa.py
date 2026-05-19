"""Stage 09 -- Magnitude-Varying Sensitivity Analysis (MV-SA).

Adapts Hadjimichael et al. (2020) MV-SA to the MOEA-FIND diagnostic
question: at each percentile of an operational hazard outcome (the
*magnitude axis*, e.g. NYC minimum reservoir storage), which
drought-hazard characteristics drive the system's response?

The factor space is the optimized MOEA-FIND objective axes (drought
characteristics from the upstream archive, same as Stage 08). The
magnitude axis is a single column from the Stage-06 metric bank.
A control uniform-random factor is appended so the magnitude-varying
noise floor is empirically visible at every percentile.

**MPI parallelization.** Percentile slices are independent — each tau
runs the full SA pipeline on the same (X, M) inputs but with a
different per-slice indicator (or window). This driver distributes
the percentile sweep across MPI ranks: rank 0 loads inputs, broadcasts
to all ranks, every rank computes a round-robin subset of percentiles,
and rank 0 gathers and writes. Run with mpirun (see
``workflows/09_magnitude_varying_sa/slurm/run_mv_sa.slurm``); a
single-rank invocation (``python run_mv_sa.py ...``) also works and
trivially reduces to serial.

Compute-only driver: writes parquet + JSON under
``outputs/09_magnitude_varying_sa/run_mv_sa/<slug>/``. Figures are
emitted by ``workflows/99_supporting_info_figures/run_mv_sa.py``.

Outputs (``outputs/09_magnitude_varying_sa/run_mv_sa/<slug>/``):

    config.json
    results/
        mv_sa_<method>.parquet      (long form: percentile x factor x cols)
        mv_sa_combined.parquet      (long form: all methods)

Usage:
    mpirun -np 19 python workflows/09_magnitude_varying_sa/run_mv_sa.py \\
        --bank  outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.csv \\
        --chars outputs/04_moea_find_single_site/run_moea_find/<src_slug>/results.json \\
        --magnitude-axis nyc_min_storage_frac \\
        --config workflows/09_magnitude_varying_sa/configs/delta_only.yaml
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.metrics.drought_metrics import PRESETS, REGISTRY, resolve_metric_set, metric_names  # noqa: E402
from src.io_paths.io import load_metric_bank, load_pareto_chars, save_experiment_config  # noqa: E402
from src.sensitivity.magnitude_varying_sa import compute_mv_sa  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402
from src.sensitivity.sensitivity import METHODS, drop_low_variance_factors  # noqa: E402
from src.io_paths.slugs import build_slug  # noqa: E402

STAGE = "09_magnitude_varying_sa"
DRIVER = "run_mv_sa"

DEFAULT_METHODS = ("delta",)
DEFAULT_MAGNITUDE_AXIS = "nyc_min_storage_frac"
DEFAULT_PERCENTILES = tuple(round(float(t), 3)
                            for t in np.linspace(0.05, 0.95, 19))


_YAML_TO_CLI = {
    "methods": "methods",
    "magnitude_axis": "magnitude_axis",
    "metric_set": "metric_set",
    "percentiles": "percentiles",
    "n_bootstrap": "n_bootstrap",
    "subsample_n": "subsample_n",
    "trace_series": "trace_series",
    "seed": "seed",
    "cv_drop_threshold": "cv_drop_threshold",
    "include_control": "include_control",
    "response_form": "response_form",
    "secondary": "secondary",
    "window_frac": "window_frac",
    "delta_num_resamples": "delta_num_resamples",
    "pawn_S": "pawn_S",
    "rbd_fast_M": "rbd_fast_M",
    "rbd_fast_num_resamples": "rbd_fast_num_resamples",
}


def _load_yaml_config(config_path: Path) -> dict:
    """Load a YAML preset into argparse-compatible defaults."""
    import yaml

    raw = yaml.safe_load(Path(config_path).read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{config_path}: top-level YAML must be a mapping")
    overrides: dict = {}
    for k, v in raw.items():
        if k in _YAML_TO_CLI:
            overrides[_YAML_TO_CLI[k]] = v
        else:
            print(f"[09] WARN: unknown YAML key {k!r} in {config_path}; ignored")
    return overrides


def _resolve_factor_set(
    chars_payload: dict,
    metric_set_arg: Optional[str],
) -> Tuple[List[str], str]:
    """Resolve the SA factor set + slug-friendly metric-set tag.

    Mirrors :func:`workflows/08_nyc_sensitivity/run_sa.py::_resolve_factor_set`.
    Default factor set comes from the upstream archive's
    ``objective_keys`` (drought-hazard characteristics MOEA-FIND
    optimized over); ``--metric-set`` overrides with a preset name or
    a comma-separated list.
    """
    if metric_set_arg is None:
        obj_keys = list(chars_payload.get("objective_keys") or [])
        if not obj_keys:
            raise SystemExit(
                "[09] no --metric-set supplied and chars JSON has no "
                "'objective_keys'. Supply --metric-set explicitly."
            )
        dropped = [k for k in obj_keys if k not in REGISTRY]
        if dropped:
            print(f"[09] dropping non-registry objective_keys entries: {dropped}")
        factor_names = [k for k in obj_keys if k in REGISTRY]
        if not factor_names:
            raise SystemExit(
                f"[09] objective_keys {obj_keys} contains no entries in "
                "drought_metrics.REGISTRY."
            )
        ms_field = chars_payload.get("metric_set")
        if isinstance(ms_field, str) and ms_field in PRESETS:
            tag = ms_field
        else:
            tag = "h" + hashlib.sha1(
                ",".join(factor_names).encode("utf-8")
            ).hexdigest()[:6]
        return factor_names, tag

    if metric_set_arg in PRESETS:
        metric_set = resolve_metric_set(metric_set_arg)
        return list(metric_names(metric_set)), metric_set_arg
    names = [n.strip() for n in metric_set_arg.split(",") if n.strip()]
    metric_set = resolve_metric_set(names)
    factor_names = list(metric_names(metric_set))
    tag = "h" + hashlib.sha1(
        ",".join(factor_names).encode("utf-8")
    ).hexdigest()[:6]
    return factor_names, tag


def _align_xym(
    chars_df: pd.DataFrame,
    bank_df: pd.DataFrame,
    factor_names: Sequence[str],
    magnitude_axis: str,
    secondary: Optional[str],
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], List[str]]:
    """Inner-join chars and bank by realization id; build (X, M, secondary, ids)."""
    common = chars_df.index.intersection(bank_df.index)
    if not len(common):
        raise SystemExit(
            f"[09] no realization ids in common between chars "
            f"({len(chars_df)}) and metric bank ({len(bank_df)})."
        )
    X = chars_df.loc[common, list(factor_names)].astype(float)
    M = pd.to_numeric(bank_df.loc[common, magnitude_axis], errors="coerce")
    if secondary is not None:
        S = pd.to_numeric(bank_df.loc[common, secondary], errors="coerce")
        finite = X.notna().all(axis=1) & M.notna() & S.notna()
        return (
            X.loc[finite].values,
            M.loc[finite].values.astype(float),
            S.loc[finite].values.astype(float),
            list(common[finite]),
        )
    finite = X.notna().all(axis=1) & M.notna()
    return (
        X.loc[finite].values,
        M.loc[finite].values.astype(float),
        None,
        list(common[finite]),
    )


def _method_kwargs(method_name: str, args) -> dict:
    """Per-method extra kwargs harvested from CLI / YAML."""
    if method_name == "delta":
        return {
            "num_resamples": args.delta_num_resamples,
            "conf_level": 0.95,
        }
    if method_name == "pawn":
        return {
            "S": args.pawn_S,
            "n_bootstrap": 0,
        }
    if method_name == "rbd_fast":
        return {
            "M": args.rbd_fast_M,
            "num_resamples": args.rbd_fast_num_resamples,
            "conf_level": 0.95,
        }
    return {}


def _percentile_grid(arg: Optional[Sequence]) -> List[float]:
    """Resolve the percentile grid from the CLI/YAML payload."""
    if arg is None or len(arg) == 0:
        return list(DEFAULT_PERCENTILES)
    return [float(t) for t in arg]


def _build_argparser(yaml_defaults: dict) -> argparse.ArgumentParser:
    """Construct the argparser; same surface for all MPI ranks."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None)

    p = argparse.ArgumentParser(
        parents=[pre],
        description="Stage 09 -- magnitude-varying sensitivity analysis (MPI).",
    )
    p.add_argument("--bank", type=Path, required=True,
                   help="Stage-06 metric_bank (parquet or csv).")
    p.add_argument("--chars", type=Path, required=True,
                   help="Upstream MOEA-FIND results.json with "
                        "'pareto_chars' and 'objective_keys'.")
    p.add_argument("--magnitude-axis", default=DEFAULT_MAGNITUDE_AXIS,
                   help="Metric-bank column whose percentiles index the "
                        "MV-SA sweep. Default: nyc_min_storage_frac.")
    p.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS),
                   choices=sorted(METHODS.keys()))
    p.add_argument("--metric-set", default=None,
                   help="Preset name or comma-list; default reads from "
                        "upstream chars JSON's objective_keys.")
    p.add_argument("--percentiles", nargs="+", type=float, default=None,
                   help=f"Percentile grid in (0, 1). Default: "
                        f"{list(DEFAULT_PERCENTILES)}")
    p.add_argument("--n-bootstrap", type=int, default=50,
                   help="Bootstrap replicates per percentile per method.")
    p.add_argument("--subsample-n", type=int, default=0,
                   help="If >0, randomly subsample to this many "
                        "realizations on rank 0 before broadcast. "
                        "0 = use all aligned realizations. Useful for "
                        "fast iteration while methodology is being "
                        "finalized; SALib delta scales O(N^2) per "
                        "internal Plischke resample so subsampling "
                        "from 3308 to 500 buys ~50x speedup.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--cv-drop-threshold", type=float, default=0.05,
                   help="Drop factors whose |std/mean| is below this CV.")
    p.add_argument("--include-control",
                   dest="include_control",
                   action="store_true",
                   default=True,
                   help="Append a uniform-random control factor "
                        "(default: enabled).")
    p.add_argument("--no-include-control",
                   dest="include_control", action="store_false",
                   help="Disable the control factor.")
    p.add_argument("--response-form",
                   choices=("exceedance", "conditional",
                            "within_trace_percentile"),
                   default="within_trace_percentile",
                   help="within_trace_percentile (default, production) = "
                        "Y_i(tau) is the tau-th percentile of trace i's "
                        "own annual series (Hadjimichael Variant 2; "
                        "requires --trace-series). exceedance = SA on "
                        "I{M <= M(tau)} (binary; degenerate at small N -- "
                        "diagnostic only). conditional = SA on "
                        "--secondary within a window of width "
                        "--window-frac centered on tau.")
    p.add_argument("--trace-series", type=Path, default=None,
                   help="Parquet of per-trace annual operational outcome "
                        "(rows = years, columns = realization ids). "
                        "Required when --response-form is "
                        "within_trace_percentile. Built by "
                        "src.hydrology.precompute_trace_series.")
    p.add_argument("--secondary", default=None,
                   help="Metric-bank column for conditional response.")
    p.add_argument("--window-frac", type=float, default=0.30,
                   help="Window fraction for conditional response.")
    p.add_argument("--delta-num-resamples", type=int, default=100)
    p.add_argument("--pawn-S", type=int, default=10)
    p.add_argument("--rbd-fast-M", type=int, default=10)
    p.add_argument("--rbd-fast-num-resamples", type=int, default=100)
    if yaml_defaults:
        p.set_defaults(**yaml_defaults)
    return p


def _round_robin(items: List, rank: int, size: int) -> List:
    """Round-robin partition of ``items`` to one MPI rank.

    Ensures balanced load when len(items) is not a multiple of size.
    """
    return [v for i, v in enumerate(items) if i % size == rank]


def main():
    # MPI bootstrap: imports are cheap if mpi4py is installed; if it
    # isn't, we want a graceful "single rank" fallback for development.
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    except ImportError:
        comm = None
        rank, size = 0, 1

    # All ranks parse args identically (cheap, deterministic). Only rank
    # 0 reads the YAML config to keep the disk footprint small; the
    # parsed defaults are broadcast.
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None)
    pre_args, _ = pre.parse_known_args()

    # All ranks load the YAML independently. This avoids a pickle-based
    # mpi4py.bcast that has been observed to fail with MPI_ERR_OTHER on
    # the Hopper MPI installation (jobs 217980-217982 on c0001 and c0006,
    # 2026-04-30) even though identical code worked on earlier jobs. The
    # YAML file is on shared storage and ~1 KB, so 19 simultaneous reads
    # are negligibly cheap and removes the bcast as a failure point.
    yaml_defaults = (
        _load_yaml_config(pre_args.config) if pre_args.config else {}
    )
    if rank != 0:
        # Suppress the per-rank "WARN: unknown YAML key" spam that
        # _load_yaml_config emits for unknown keys; only rank 0 reports.
        pass

    p = _build_argparser(yaml_defaults or {})
    args = p.parse_args()

    if args.response_form == "conditional" and not args.secondary:
        if rank == 0:
            print("[09] --response-form conditional requires --secondary "
                  "<metric-bank-column>.", flush=True)
        raise SystemExit(2)
    if (args.response_form == "within_trace_percentile"
            and not args.trace_series):
        if rank == 0:
            print("[09] --response-form within_trace_percentile requires "
                  "--trace-series <path-to-parquet>.", flush=True)
        raise SystemExit(2)

    # ------------------------------------------------------------------
    # Every rank loads + aligns inputs independently.
    #
    # We deliberately do not use mpi4py.bcast to broadcast the payload
    # from rank 0: the same code path failed with MPI_ERR_OTHER on jobs
    # 217980-217982 (multiple Hopper nodes, 2026-04-30) even though
    # identical bcast calls worked on earlier jobs. Per-rank loads are
    # cheap (~5 s of shared-FS I/O for chars JSON + metric bank CSV +
    # trace-series parquet) and remove pickle-bcast as a failure point.
    # Only rank 0 writes config.json to avoid race conditions.
    # ------------------------------------------------------------------
    if rank == 0:
        print(f"[09 r{rank}/{size}] loading chars from {args.chars}",
              flush=True)
    chars_payload = json.loads(Path(args.chars).read_text())
    chars_df = load_pareto_chars(args.chars)

    if rank == 0:
        print(f"[09 r{rank}/{size}] loading metric bank from {args.bank}",
              flush=True)
    bank_df = load_metric_bank(args.bank)

    factor_names, metric_set_tag = _resolve_factor_set(
        chars_payload, args.metric_set,
    )
    if rank == 0:
        print(f"[09 r{rank}/{size}] factor set ({metric_set_tag}): "
              f"{factor_names}", flush=True)

    missing = [f for f in factor_names if f not in chars_df.columns]
    if missing:
        raise SystemExit(
            f"[09] factor names {missing} are not columns of "
            f"pareto_chars. Available: {list(chars_df.columns)}"
        )
    if args.magnitude_axis not in bank_df.columns:
        raise SystemExit(
            f"[09] magnitude axis {args.magnitude_axis!r} not in "
            f"metric bank. Available: "
            f"{sorted(bank_df.columns.tolist())}"
        )
    if (args.secondary is not None
            and args.secondary not in bank_df.columns):
        raise SystemExit(
            f"[09] secondary {args.secondary!r} not in metric bank. "
            f"Available: {sorted(bank_df.columns.tolist())}"
        )

    factor_view = chars_df[list(factor_names)]
    factor_view_kept, dropped_factors = drop_low_variance_factors(
        factor_view, cv_threshold=args.cv_drop_threshold,
    )
    if rank == 0 and dropped_factors:
        print(f"[09 r{rank}/{size}] dropping low-variance factors: "
              f"{dropped_factors}", flush=True)
    factor_names = list(factor_view_kept.columns)
    if not factor_names:
        raise SystemExit("[09] all factors dropped; nothing to analyze.")

    X_arr, M_arr, S_arr, kept_ids = _align_xym(
        chars_df, bank_df, factor_names, args.magnitude_axis,
        secondary=args.secondary,
    )

    # Load + align per-trace series for within_trace_percentile mode.
    trace_arr = None
    trace_n_years = None
    if args.trace_series:
        if rank == 0:
            print(f"[09 r{rank}/{size}] loading trace series from "
                  f"{args.trace_series}", flush=True)
        trace_df = pd.read_parquet(args.trace_series)
        # Columns are realization ids, rows are years; orient so each
        # row of trace_arr corresponds to one realization.
        trace_df.columns = trace_df.columns.astype(str)
        ts_aligned = trace_df.reindex(columns=kept_ids).T
        missing_ts = ts_aligned.isna().any(axis=1)
        if missing_ts.any():
            if rank == 0:
                print(f"[09 r{rank}/{size}] WARN: {int(missing_ts.sum())} "
                      f"realization ids absent from trace_series; "
                      f"dropping.", flush=True)
            keep_mask = ~missing_ts.values
            X_arr = X_arr[keep_mask]
            M_arr = M_arr[keep_mask]
            if S_arr is not None:
                S_arr = S_arr[keep_mask]
            kept_ids = [kept_ids[i] for i in range(len(kept_ids))
                        if keep_mask[i]]
            ts_aligned = ts_aligned.loc[~missing_ts]
        trace_arr = ts_aligned.values.astype(float)
        trace_n_years = int(trace_arr.shape[1])
        if rank == 0:
            print(f"[09 r{rank}/{size}] trace_series aligned: "
                  f"{trace_arr.shape[0]} traces x "
                  f"{trace_n_years} year-bins per trace",
                  flush=True)

    n_aligned = int(X_arr.shape[0])
    if args.subsample_n and args.subsample_n < n_aligned:
        sub_rng = np.random.default_rng(args.seed)
        sub_idx = sub_rng.choice(
            n_aligned, size=args.subsample_n, replace=False,
        )
        sub_idx.sort()
        X_arr = X_arr[sub_idx]
        M_arr = M_arr[sub_idx]
        if S_arr is not None:
            S_arr = S_arr[sub_idx]
        if trace_arr is not None:
            trace_arr = trace_arr[sub_idx]
        kept_ids = [kept_ids[i] for i in sub_idx]
        if rank == 0:
            print(f"[09 r{rank}/{size}] subsampled "
                  f"{args.subsample_n}/{n_aligned} realizations "
                  f"(seed={args.seed})", flush=True)
    n_realizations = int(X_arr.shape[0])
    if rank == 0:
        print(f"[09 r{rank}/{size}] aligned n_realizations={n_realizations}, "
              f"M range=[{M_arr.min():.4g}, {M_arr.max():.4g}]",
              flush=True)

    percentiles = _percentile_grid(args.percentiles)

    upstream_slug = (chars_payload.get("variant")
                     or Path(args.chars).parent.name)
    slug = build_slug(
        "mvsa",
        src=upstream_slug,
        axis=args.magnitude_axis,
        resp=args.response_form,
        metric_set=metric_set_tag,
        methods="-".join(sorted(args.methods)),
        n_perc=len(percentiles),
        n_sub=(args.subsample_n if args.subsample_n else None),
        s=args.seed,
    )
    out = stage_output_dir(STAGE, DRIVER, slug)
    results_dir = out / "results"
    if rank == 0:
        results_dir.mkdir(parents=True, exist_ok=True)
        print(f"[09 r{rank}/{size}] variant: {slug}", flush=True)
        print(f"[09 r{rank}/{size}] output:  {out}", flush=True)
        save_experiment_config(out, {
            "script": "workflows/09_magnitude_varying_sa/run_mv_sa.py",
            "stage": STAGE,
            "driver": DRIVER,
            "variant": slug,
            "bank": str(args.bank),
            "chars": str(args.chars),
            "methods": list(args.methods),
            "magnitude_axis": args.magnitude_axis,
            "metric_set_tag": metric_set_tag,
            "factor_names": list(factor_names),
            "dropped_factors": dropped_factors,
            "percentiles": list(percentiles),
            "n_bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "cv_drop_threshold": args.cv_drop_threshold,
            "include_control": bool(args.include_control),
            "response_form": args.response_form,
            "secondary": args.secondary,
            "window_frac": float(args.window_frac),
            "n_realizations": n_realizations,
            "n_aligned": n_aligned,
            "subsample_n": int(args.subsample_n),
            "trace_series_path": (str(args.trace_series)
                                  if args.trace_series else None),
            "trace_n_years": trace_n_years,
            "mpi_size": size,
        })

    # ------------------------------------------------------------------
    # Each rank computes its round-robin subset of percentiles, for
    # every requested method, then writes its partial parquet to disk.
    # Rank 0 polls the filesystem for completion sentinels and
    # concatenates without using mpi4py.gather (which is also pickle-
    # bcast-based and would fail under the same MPI condition that
    # broke yaml_defaults bcast on jobs 217980-217982).
    # ------------------------------------------------------------------
    my_percentiles = _round_robin(percentiles, rank, size)
    print(f"[09 r{rank}/{size}] my_percentiles ({len(my_percentiles)}/"
          f"{len(percentiles)}): {my_percentiles}", flush=True)

    # Per-rank scratch dir for partial parquets + sentinel files.
    # All ranks mkdir(exist_ok=True) -- avoids a race where non-rank-0
    # ranks try to write a sentinel before rank 0 has created the dir.
    partial_dir = results_dir / "_partials"
    partial_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    for m_idx, method_name in enumerate(args.methods, 1):
        kw = _method_kwargs(method_name, args)
        if rank == 0:
            print(f"[09 r0/{size}] method {m_idx}/{len(args.methods)} = "
                  f"{method_name} (elapsed={time.time()-t0:.1f}s)",
                  flush=True)

        if my_percentiles:
            t_local = time.time()
            local_df = compute_mv_sa(
                X_arr, M_arr, factor_names,
                method=method_name,
                response_form=args.response_form,
                percentiles=my_percentiles,
                secondary=S_arr,
                trace_series=trace_arr,
                window_frac=args.window_frac,
                n_bootstrap=args.n_bootstrap,
                seed=args.seed,
                include_control=args.include_control,
                method_kwargs=kw,
            )
            local_df["axis_column"] = args.magnitude_axis
            print(f"[09 r{rank}/{size}]   {method_name}: "
                  f"{len(my_percentiles)} percentiles in "
                  f"{time.time()-t_local:.1f}s", flush=True)
        else:
            local_df = pd.DataFrame()

        # Each rank writes its partial parquet + a sentinel.
        partial_path = partial_dir / f"mv_sa_{method_name}_r{rank:02d}.parquet"
        sentinel_path = partial_dir / f"mv_sa_{method_name}_r{rank:02d}.done"
        if not local_df.empty:
            local_df.to_parquet(partial_path)
        sentinel_path.touch()

        # Rank 0 waits for all ranks' sentinels, then concatenates.
        if rank == 0:
            expected = {
                partial_dir / f"mv_sa_{method_name}_r{r:02d}.done"
                for r in range(size)
            }
            t_wait = time.time()
            while True:
                missing = [p for p in expected if not p.exists()]
                if not missing:
                    break
                if time.time() - t_wait > 1800:  # 30-min hard wait cap
                    print(f"[09 r0/{size}] WARN: timed out waiting for "
                          f"{len(missing)} ranks; concatenating partials.",
                          flush=True)
                    break
                time.sleep(1.0)
            partials = sorted(
                partial_dir.glob(f"mv_sa_{method_name}_r*.parquet")
            )
            non_empty = [pd.read_parquet(p) for p in partials]
            df_m = (pd.concat(non_empty, ignore_index=True)
                    if non_empty else pd.DataFrame())
            df_m = (df_m.sort_values(["percentile", "factor"])
                        .reset_index(drop=True))
            out_path = results_dir / f"mv_sa_{method_name}.parquet"
            df_m.to_parquet(out_path)
            print(f"[09 r0/{size}]   wrote {out_path.name} "
                  f"({len(df_m)} rows)", flush=True)

    if rank == 0:
        # Combined parquet across methods (re-read cheaply from disk).
        combined_frames = [
            pd.read_parquet(results_dir / f"mv_sa_{m}.parquet")
            for m in args.methods
            if (results_dir / f"mv_sa_{m}.parquet").exists()
        ]
        if combined_frames:
            combined = pd.concat(combined_frames, ignore_index=True)
            combined.to_parquet(results_dir / "mv_sa_combined.parquet")
        print(f"[09 r0/{size}] computed all methods in "
              f"{time.time() - t0:.1f}s", flush=True)
        print(f"[09 r0/{size}] done. results -> {results_dir}", flush=True)


if __name__ == "__main__":
    main()
