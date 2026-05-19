"""Kirsch-wrapper mapping geometry diagnostic (SI C, compute).

Sweeps a single DV coordinate over [0, 1] for both wrapper modes (index
vs residual), with all other DVs held at 0.5, and records the resulting
flow response. Plus a synthetic per-representative-DV seasonal profile
used by the smoothness figure.

Compute only -- writes numerical artifacts under
``outputs/02_calibration/wrapper_geometry/``:
    - config.json
    - geometry_summary.json
    - sweep.npz   (per-mode p_grid, flow_at_k, annual_totals, k, month_idx,
                   plus rep_profiles for representative DV values)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.hydrology.kirsch_utils import build_kirsch_generator  # noqa: E402
from src.hydrology.kirsch_wrapper import KirschBorgWrapper  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "02_calibration"
DRIVER = "wrapper_geometry"

SWEEP_K_FRACTION = 0.5
REP_DVS = (0.1, 0.5, 0.9)


def run_sweep(wrapper, k: int, n_sweep: int) -> dict:
    n_dvs = wrapper.n_dvs
    n_years = wrapper.n_years_out
    v0 = np.full(n_dvs, 0.5)
    month_idx = k % 12

    p_grid = np.linspace(0.0, 1.0, n_sweep)
    flow_at_k = np.zeros(n_sweep)
    annual_totals = np.zeros((n_sweep, n_years))

    for i, p in enumerate(p_grid):
        v = v0.copy()
        v[k] = p
        synthetic = wrapper.generate(v)
        water_col = (month_idx + 3) % 12
        flow_at_k[i] = synthetic[:, water_col].mean()
        annual_totals[i, :] = synthetic.sum(axis=1)

    distinct_vals = len(np.unique(np.round(flow_at_k, decimals=6)))

    # Representative DV seasonal profiles (median-year shape per rep DV).
    rep_profiles = np.zeros((len(REP_DVS), 12))
    for s_idx, p_rep in enumerate(REP_DVS):
        v = v0.copy()
        v[k] = p_rep
        syn = wrapper.generate(v)
        year_totals = syn.sum(axis=1)
        med_yr = int(np.argsort(year_totals)[n_years // 2])
        rep_profiles[s_idx] = syn[med_yr, :]

    return {
        "p_grid": p_grid,
        "flow_at_k": flow_at_k,
        "annual_totals": annual_totals,
        "annual_mean": annual_totals.mean(axis=1),
        "distinct_count": distinct_vals,
        "k": k,
        "month_idx": month_idx,
        "rep_profiles": rep_profiles,
    }


def _flow_change_stats(p_grid: np.ndarray, flow: np.ndarray) -> dict:
    dp = np.diff(p_grid)
    df = np.abs(np.diff(flow))
    rate = df / np.where(dp > 0, dp, np.nan)
    return {
        "mean_change_per_unit_dv": float(np.nanmean(rate)),
        "max_change_per_unit_dv": float(np.nanmax(rate)),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-sweep-points", type=int, default=400)
    p.add_argument("--n-years", type=int, default=20)
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    (out_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "seed": args.seed,
        "n_sweep_points": args.n_sweep_points,
        "n_years": args.n_years,
        "rep_dvs": list(REP_DVS),
    }, indent=2))

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)
    n_years_hist = monthly_2d.shape[0]
    print(f"[02/wrapper_geometry] n_years_hist={n_years_hist}")

    gen = build_kirsch_generator(monthly_2d)
    wrapper_index = KirschBorgWrapper(gen, mode="index", n_years_out=args.n_years)
    wrapper_residual = KirschBorgWrapper(gen, mode="residual", n_years_out=args.n_years)

    k_residual = int(SWEEP_K_FRACTION * wrapper_residual.n_dvs)
    k_index = min(k_residual, wrapper_index.n_dvs - 1)
    print(f"  sweeping residual k={k_residual} (month {k_residual % 12}), "
          f"index k={k_index} (month {k_index % 12})")

    print("  running index sweep ...")
    res_index = run_sweep(wrapper_index, k_index, args.n_sweep_points)
    print("  running residual sweep ...")
    res_residual = run_sweep(wrapper_residual, k_residual, args.n_sweep_points)

    np.savez(
        out_dir / "sweep.npz",
        index_p_grid=res_index["p_grid"],
        index_flow_at_k=res_index["flow_at_k"],
        index_annual_totals=res_index["annual_totals"],
        index_annual_mean=res_index["annual_mean"],
        index_rep_profiles=res_index["rep_profiles"],
        index_k=res_index["k"],
        index_month_idx=res_index["month_idx"],
        index_distinct_count=res_index["distinct_count"],
        residual_p_grid=res_residual["p_grid"],
        residual_flow_at_k=res_residual["flow_at_k"],
        residual_annual_totals=res_residual["annual_totals"],
        residual_annual_mean=res_residual["annual_mean"],
        residual_rep_profiles=res_residual["rep_profiles"],
        residual_k=res_residual["k"],
        residual_month_idx=res_residual["month_idx"],
        residual_distinct_count=res_residual["distinct_count"],
        n_years_hist=n_years_hist,
        rep_dvs=np.array(REP_DVS),
    )
    print(f"  wrote {out_dir / 'sweep.npz'}")

    summary = {"n_years_hist": n_years_hist, "modes": {}}
    for mode, res in (("index", res_index), ("residual", res_residual)):
        stats = _flow_change_stats(res["p_grid"], res["flow_at_k"])
        summary["modes"][mode] = {
            "distinct_flow_values": int(res["distinct_count"]),
            "swept_coordinate_k": int(res["k"]),
            "controlled_month_idx": int(res["month_idx"]),
            **stats,
        }
    (out_dir / "geometry_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {out_dir / 'geometry_summary.json'}")


if __name__ == "__main__":
    main()
