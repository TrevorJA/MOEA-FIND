"""Kirsch-wrapper fidelity diagnostic (SI C, compute).

Compares four ensembles of T-year traces at the trace level:

    historical        -- overlapping T-year blocks of the historical record.
    kirsch_baseline   -- traces from SynHydro's KirschGenerator.generate()
                         (no DV wrapper). The "natural" Kirsch distribution.
    index             -- KirschBorgWrapper(mode="index") under U[0,1] DVs.
    residual          -- KirschBorgWrapper(mode="residual") under U[0,1] DVs.

Phase A asks whether each wrapper reproduces the underlying Kirsch
distribution when DVs are drawn uniformly. Phase B converts each trace
into drought characteristics and overlays mode clouds on the historical
T-block cloud, an in-process kirsch_baseline cloud, and the stage-03
Kirsch library (when present).

Compute only -- writes ensembles and summary JSON under
``outputs/02_calibration/wrapper_fidelity/`` for the paired plotting
driver to consume.

Outputs:
    - config.json
    - ensembles_2d.npz          (historical, kirsch_baseline, index, residual,
                                 each (n_traces, n_years, 12))
    - dvs.npz                   (index, residual U[0,1] DV matrices)
    - drought_clouds.npz        (mode clouds, hist_cloud, baseline_cloud,
                                 exp05_cloud, objective_keys, objective_labels)
    - fidelity_summary.json
    - coverage_summary.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment_utils import (  # noqa: E402
    compute_historical_ssi_chars,
    prepare_data,
)
from src.historical_blocks import (  # noqa: E402
    compute_historical_block_chars,
    resample_historical_blocks,
    resample_historical_blocks_2d,
)
from src.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.objectives import (  # noqa: E402
    compute_ssi_drought_characteristics,
    flows_to_series,
)
from src.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "wrapper_fidelity"

ENS_ORDER: Tuple[str, ...] = ("historical", "kirsch_baseline", "index", "residual")


# ----------------------------------------------------------------------------
# Generator helpers
# ----------------------------------------------------------------------------

def _baseline_traces(
    kirsch_gen,
    n_traces: int,
    n_years: int,
    seed: int,
) -> np.ndarray:
    """Generate kirsch_baseline traces via SynHydro's native generate().

    Returns array of shape (n_traces, n_years, 12) in water-year order.
    """
    ensemble = kirsch_gen.generate(
        n_realizations=n_traces,
        n_years=n_years,
        seed=seed,
    )
    out = np.empty((n_traces, n_years, 12), dtype=float)
    for i, rid in enumerate(sorted(ensemble.data_by_realization.keys())):
        trace_df = ensemble.data_by_realization[rid]
        trace_1d = trace_df.values.flatten()
        n_yrs = len(trace_1d) // 12
        trace_2d = trace_1d[: n_yrs * 12].reshape(n_yrs, 12)
        # SynHydro emits calendar-year order; roll to water-year (Oct first).
        out[i] = np.roll(trace_2d, 3, axis=1)
    return out


def _wrapper_chunk_worker(
    mode: str,
    n_traces_chunk: int,
    n_years: int,
    seed_chunk: int,
    monthly_2d: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    kirsch_gen = build_kirsch_generator(monthly_2d)
    wrapper = KirschBorgWrapper(kirsch_gen, mode=mode, n_years_out=n_years)
    rng = np.random.default_rng(seed_chunk)
    dvs = rng.uniform(0.0, 1.0, size=(n_traces_chunk, wrapper.n_dvs))
    traces = np.empty((n_traces_chunk, n_years, 12), dtype=float)
    for i in range(n_traces_chunk):
        traces[i] = np.asarray(wrapper.generate(dvs[i]), dtype=float)
    return traces, dvs


def _wrapper_traces(
    kirsch_gen,
    mode: str,
    n_traces: int,
    n_years: int,
    seed: int,
    monthly_2d: np.ndarray,
    workers: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Returns (traces array (n_traces, n_years, 12), dvs (n_traces, n_dvs))."""
    if workers <= 1:
        wrapper = KirschBorgWrapper(kirsch_gen, mode=mode, n_years_out=n_years)
        rng = np.random.default_rng(seed)
        dvs = rng.uniform(0.0, 1.0, size=(n_traces, wrapper.n_dvs))
        traces = np.empty((n_traces, n_years, 12), dtype=float)
        for i in range(n_traces):
            traces[i] = np.asarray(wrapper.generate(dvs[i]), dtype=float)
        return traces, dvs

    chunk_sizes = [n_traces // workers] * workers
    for i in range(n_traces % workers):
        chunk_sizes[i] += 1
    seeds = [seed + 1_000_003 * i for i in range(workers)]

    traces_chunks: List[np.ndarray] = [None] * workers  # type: ignore
    dvs_chunks: List[np.ndarray] = [None] * workers  # type: ignore
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                _wrapper_chunk_worker, mode, chunk_sizes[i], n_years,
                seeds[i], monthly_2d,
            ): i
            for i in range(workers) if chunk_sizes[i] > 0
        }
        for fut in as_completed(futures):
            i = futures[fut]
            tr, dv = fut.result()
            traces_chunks[i] = tr
            dvs_chunks[i] = dv

    traces_all = np.concatenate([c for c in traces_chunks if c is not None], axis=0)
    dvs_all = np.concatenate([c for c in dvs_chunks if c is not None], axis=0)
    return traces_all, dvs_all


# ----------------------------------------------------------------------------
# Statistical primitives (used by summary JSON only; plot driver re-derives)
# ----------------------------------------------------------------------------

def _acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    n = len(x)
    mean = np.mean(x)
    var = np.var(x)
    if var < 1e-12:
        return np.zeros(max_lag + 1)
    out = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        out[lag] = np.mean((x[: n - lag] - mean) * (x[lag:] - mean)) / var
    return out


def _cross_month_corr(traces_2d: np.ndarray) -> np.ndarray:
    mats = []
    for t in traces_2d:
        if t.shape[0] < 2:
            continue
        mats.append(np.corrcoef(t, rowvar=False))
    return np.mean(np.array(mats), axis=0) if mats else np.zeros((12, 12))


def _annual_totals(trace_2d: np.ndarray) -> np.ndarray:
    return trace_2d.sum(axis=1)


def _annual_min_monthly(trace_2d: np.ndarray) -> np.ndarray:
    return trace_2d.min(axis=1)


def _ks_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import ks_2samp
    return float(ks_2samp(a, b).statistic)


def _ad_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import anderson_ksamp
    try:
        return float(anderson_ksamp([a, b]).statistic)
    except Exception:
        return float("nan")


# ----------------------------------------------------------------------------
# Drought-space (Phase B)
# ----------------------------------------------------------------------------

def _drought_chars_chunk_worker(
    traces_1d: List[np.ndarray],
    ssi_calc,
    objective_keys: Tuple[str, ...],
) -> np.ndarray:
    out = np.zeros((len(traces_1d), len(objective_keys)))
    for i, t in enumerate(traces_1d):
        series = flows_to_series(t, start_date="2100-01-01")
        ssi = ssi_calc.transform(series)
        chars = compute_ssi_drought_characteristics(ssi, monthly_flows=t)
        for j, k in enumerate(objective_keys):
            out[i, j] = float(chars.get(k, np.nan))
    return out


def _drought_chars_for_traces(
    traces_1d: List[np.ndarray],
    ssi_calc,
    objective_keys: Tuple[str, ...],
    workers: int = 1,
) -> np.ndarray:
    n = len(traces_1d)
    if workers <= 1 or n < workers * 4:
        return _drought_chars_chunk_worker(traces_1d, ssi_calc, objective_keys)

    chunk_size = n // workers
    chunks = [
        traces_1d[i * chunk_size : (i + 1) * chunk_size]
        for i in range(workers - 1)
    ]
    chunks.append(traces_1d[(workers - 1) * chunk_size :])

    parts: List[np.ndarray] = [None] * workers  # type: ignore
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                _drought_chars_chunk_worker, chunks[i], ssi_calc, objective_keys,
            ): i
            for i in range(workers) if len(chunks[i]) > 0
        }
        for fut in as_completed(futures):
            parts[futures[fut]] = fut.result()
    return np.concatenate([p for p in parts if p is not None], axis=0)


def _hull_volume(points: np.ndarray) -> float:
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return float("nan")
    pts = points[~np.isnan(points).any(axis=1)]
    if len(pts) < pts.shape[1] + 1:
        return float("nan")
    try:
        return float(ConvexHull(pts).volume)
    except Exception:
        return float("nan")


def _fraction_inside_hull(points: np.ndarray, hull_points: np.ndarray) -> float:
    try:
        from scipy.spatial import Delaunay
    except ImportError:
        return float("nan")
    pts = points[~np.isnan(points).any(axis=1)]
    hpts = hull_points[~np.isnan(hull_points).any(axis=1)]
    if len(hpts) < hpts.shape[1] + 1 or len(pts) == 0:
        return float("nan")
    try:
        tri = Delaunay(hpts)
        return float(np.mean(tri.find_simplex(pts) >= 0))
    except Exception:
        return float("nan")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-traces", type=int, default=1000,
                   help="Traces per ensemble (wrappers and kirsch_baseline).")
    p.add_argument("--n-years", type=int, default=20,
                   help="Trace length in water years.")
    p.add_argument("--ssi", type=int, default=3, choices=[1, 3, 6, 12])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--kirsch-library-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "03_kirsch_library",
                   help="Stage-03 Kirsch library directory (Phase B reference cloud).")
    p.add_argument("--skip-phase-b", action="store_true",
                   help="Skip the drought-space coverage section (Phase B).")
    p.add_argument("--metric-set", default="primary",
                   help="Drought metric set name from src.drought_metrics.PRESETS.")
    p.add_argument("--workers", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
                   help="Parallel worker count.")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/wrapper_fidelity] output_dir={out_dir}")
    (out_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "n_traces": args.n_traces, "n_years": args.n_years,
        "ssi": args.ssi, "seed": args.seed,
        "metric_set": args.metric_set,
        "kirsch_library_dir": str(args.kirsch_library_dir),
        "skip_phase_b": bool(args.skip_phase_b),
    }, indent=2))

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    print(f"  historical: {monthly_2d.shape[0]} water years")

    hist_blocks_1d = resample_historical_blocks(monthly_1d, T_years=args.n_years, stride=1)
    hist_blocks_2d = resample_historical_blocks_2d(monthly_2d, T_years=args.n_years, stride=1)
    print(f"  historical T-blocks: n={len(hist_blocks_1d)}")
    hist_2d_arr = np.stack(hist_blocks_2d, axis=0)

    print(f"  fitting Kirsch generator ...")
    kirsch_gen = build_kirsch_generator(monthly_2d)

    print(f"  generating kirsch_baseline (n={args.n_traces}) ...")
    baseline_2d = _baseline_traces(
        kirsch_gen, args.n_traces, args.n_years, args.seed,
    )

    workers = max(1, int(args.workers))
    print(f"  workers={workers}")

    print(f"  generating index-wrapper traces ...")
    index_2d, index_dvs = _wrapper_traces(
        kirsch_gen, "index", args.n_traces, args.n_years, args.seed + 1,
        monthly_2d=monthly_2d, workers=workers,
    )

    print(f"  generating residual-wrapper traces ...")
    resid_2d, resid_dvs = _wrapper_traces(
        kirsch_gen, "residual", args.n_traces, args.n_years, args.seed + 2,
        monthly_2d=monthly_2d, workers=workers,
    )

    np.savez(
        out_dir / "ensembles_2d.npz",
        historical=hist_2d_arr,
        kirsch_baseline=baseline_2d,
        index=index_2d,
        residual=resid_2d,
    )
    np.savez(out_dir / "dvs.npz", index=index_dvs, residual=resid_dvs)
    print(f"  wrote ensembles_2d.npz, dvs.npz")

    # ------------------------------------------------------------------
    # Fidelity summary (per-ensemble scalar diagnostics)
    # ------------------------------------------------------------------
    ensembles_2d_arr: Dict[str, np.ndarray] = {
        "historical": hist_2d_arr,
        "kirsch_baseline": baseline_2d,
        "index": index_2d,
        "residual": resid_2d,
    }
    hist_flat_by_month = [hist_2d_arr[:, :, m].flatten() for m in range(12)]
    summary: Dict = {
        "config": {
            "n_traces": args.n_traces, "n_years": args.n_years,
            "seed": args.seed,
        },
        "n_hist_blocks": len(hist_blocks_1d),
    }
    summary_by_ens: Dict[str, Dict] = {}
    cxm_hist = _cross_month_corr(hist_2d_arr)
    for ens in ENS_ORDER:
        traces_2d = ensembles_2d_arr[ens]
        if traces_2d.shape[0] == 0:
            continue
        traces_1d = traces_2d.reshape(traces_2d.shape[0], -1)
        annual_totals = traces_2d.sum(axis=2).flatten()
        annual_mins = traces_2d.min(axis=2).flatten()
        acfs = np.array([_acf(t, 12) for t in traces_1d])
        cxm = _cross_month_corr(traces_2d)
        flat_by_month = [traces_2d[:, :, m].flatten() for m in range(12)]
        per_month_ks = [_ks_two_sample(flat_by_month[m], hist_flat_by_month[m])
                        for m in range(12)]
        per_month_ad = [_ad_two_sample(flat_by_month[m], hist_flat_by_month[m])
                        for m in range(12)]
        summary_by_ens[ens] = {
            "n_traces": int(traces_2d.shape[0]),
            "annual_mean": float(annual_totals.mean()),
            "annual_cv": float(annual_totals.std(ddof=1) / annual_totals.mean()),
            "annual_min_median": float(np.median(annual_mins)),
            "lag1_ac_median": float(np.median(acfs[:, 1])),
            "lag12_ac_median": float(np.median(acfs[:, 12])),
            "cross_month_frob_vs_hist": float(np.linalg.norm(cxm - cxm_hist, ord="fro")),
            "per_month_ks_vs_hist": per_month_ks,
            "per_month_ad_vs_hist": per_month_ad,
            "per_month_ks_max": float(max(per_month_ks)),
            "per_month_ad_max": float(max(per_month_ad)),
        }
    summary["ensembles"] = summary_by_ens
    (out_dir / "fidelity_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    print(f"  wrote {out_dir / 'fidelity_summary.json'}")

    # ------------------------------------------------------------------
    # Phase B -- drought-space coverage
    # ------------------------------------------------------------------
    if args.skip_phase_b:
        print(f"  Phase B skipped.")
        return

    print(f"  Phase B: drought characterisation ...")
    from src.drought_metrics import metric_names, resolve_metric_set
    metric_set = resolve_metric_set(args.metric_set)
    objective_keys = metric_names(metric_set)
    objective_labels = tuple(f"{m.label} ({m.units})" for m in metric_set)
    print(f"  metric set: {args.metric_set} -> {objective_keys}")
    _, ssi_calc, _ = compute_historical_ssi_chars(monthly_1d, args.ssi)

    index_1d = [index_2d[i].flatten() for i in range(index_2d.shape[0])]
    resid_1d = [resid_2d[i].flatten() for i in range(resid_2d.shape[0])]
    baseline_1d = [baseline_2d[i].flatten() for i in range(baseline_2d.shape[0])]

    mode_clouds: Dict[str, np.ndarray] = {
        "index": _drought_chars_for_traces(index_1d, ssi_calc, objective_keys,
                                           workers=workers),
        "residual": _drought_chars_for_traces(resid_1d, ssi_calc, objective_keys,
                                              workers=workers),
    }
    in_proc_baseline_cloud = _drought_chars_for_traces(
        baseline_1d, ssi_calc, objective_keys, workers=workers,
    )
    hist_cloud = compute_historical_block_chars(
        monthly_1d, T_years=args.n_years, ssi_calc=ssi_calc,
        objective_keys=objective_keys, stride=1,
    )

    exp05_cloud: np.ndarray = np.zeros((0, len(objective_keys)))
    exp05_path = args.kirsch_library_dir / "characteristics.npz"
    if exp05_path.exists():
        try:
            klib = np.load(exp05_path, allow_pickle=True)
            k_all_keys = [str(k) for k in klib["all_keys"]]
            k_all_vals = klib["all_values"]
            obj_indices: List[int] = []
            for ok in objective_keys:
                obj_indices.append(k_all_keys.index(ok) if ok in k_all_keys else -1)
            if all(j >= 0 for j in obj_indices):
                exp05_cloud = k_all_vals[:, obj_indices].astype(float)
                print(f"  loaded Kirsch library: n={exp05_cloud.shape[0]}")
            else:
                missing = [ok for ok, j in zip(objective_keys, obj_indices) if j < 0]
                print(f"  Kirsch library missing keys {missing}; skipping.")
        except Exception as exc:
            print(f"  WARNING loading Kirsch library: {exc}")
    else:
        print(f"  Kirsch library not found at {exp05_path}; skipping reference cloud.")

    np.savez(
        out_dir / "drought_clouds.npz",
        index=mode_clouds["index"],
        residual=mode_clouds["residual"],
        kirsch_baseline=in_proc_baseline_cloud,
        historical=hist_cloud,
        exp05=exp05_cloud,
        objective_keys=np.array(list(objective_keys)),
        objective_labels=np.array(list(objective_labels)),
    )
    print(f"  wrote {out_dir / 'drought_clouds.npz'}")

    cov: Dict = {"objective_keys": list(objective_keys)}
    per_mode_cov: Dict[str, Dict] = {}
    ref_hull_2d = exp05_cloud[:, :2] if exp05_cloud.shape[0] > 0 else hist_cloud[:, :2]
    ref_hull_3d = exp05_cloud if exp05_cloud.shape[0] > 0 else hist_cloud
    for mode in ("index", "residual"):
        cloud = mode_clouds[mode]
        if cloud.size == 0:
            continue
        per_mode_cov[mode] = {
            "n": int(cloud.shape[0]),
            "per_axis_min": [float(x) for x in np.nanmin(cloud, axis=0)],
            "per_axis_max": [float(x) for x in np.nanmax(cloud, axis=0)],
            "hull_volume_2d": _hull_volume(cloud[:, :2]),
            "hull_volume_3d": _hull_volume(cloud),
            "frac_inside_exp05_hull_2d": _fraction_inside_hull(cloud[:, :2], ref_hull_2d),
            "frac_inside_exp05_hull_3d": _fraction_inside_hull(cloud, ref_hull_3d),
        }
    cov["per_mode"] = per_mode_cov
    cov["in_proc_baseline_n"] = int(in_proc_baseline_cloud.shape[0])
    cov["exp05_library_n"] = int(exp05_cloud.shape[0])
    cov["historical_block_n"] = int(hist_cloud.shape[0])
    (out_dir / "coverage_summary.json").write_text(
        json.dumps(cov, indent=2, default=str)
    )
    print(f"  wrote {out_dir / 'coverage_summary.json'}")
    print(f"  done.")


if __name__ == "__main__":
    main()
