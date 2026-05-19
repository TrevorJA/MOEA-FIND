"""Analytic epsilon x NFE sensitivity sweep (manuscript Fig 3).

Coarse 4 eps x 3 NFE x 3 seeds = 36 cells, run as a SLURM job array on
HPC. Each cell runs MM Borg MOEA via MPI; the array index uniquely
selects (epsilon, NFE, seed). Aggregate mode pools the per-cell JSONs
into a summary table and a recommendation.

Outputs under ``outputs/01_analytic_validation/eps_nfe_sweep/``:
    - config.json
    - cells/cell_<id>.json    (one per array task)
    - aggregate.json          (after --mode aggregate)

The companion plotting driver
``workflows/99_supporting_info_figures/eps_nfe_sweep.py`` reads
``aggregate.json`` and renders the heatmap figure.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from src.discovery.analysis import coverage_metrics  # noqa: E402
from src.optimization.borg_runner import run_optimization  # noqa: E402
from src.metrics.objectives import analytic_objectives  # noqa: E402
from src.io_paths.paths import stage_output_dir  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "eps_nfe_sweep"

K = 3
DV_LO, DV_HI = -3.0, 3.0
ANTI_IDEAL = np.array([3.0, 3.0, 3.0])

# Coarsened grid: 4 eps x 3 NFE x 3 seeds = 36 cells. Each cell is one
# MM Borg MPI run; the full array fits in well under the 400-core budget.
DEFAULT_EPSILONS = [0.10, 0.15, 0.20, 0.30]
DEFAULT_NFES = [5_000, 20_000, 50_000]
DEFAULT_SEEDS = [42, 17, 91]
DEFAULT_METRIC_CAP = 0  # 0 = no L2*-discrepancy subsample (always full on HPC)


def _cell_index(eps_i, nfe_i, seed_i, n_nfe, n_seeds):
    return eps_i * n_nfe * n_seeds + nfe_i * n_seeds + seed_i


def _unpack_cell(task_id, eps_list, nfe_list, seed_list):
    n_nfe = len(nfe_list); n_seeds = len(seed_list)
    eps_i, rem = divmod(task_id, n_nfe * n_seeds)
    nfe_i, seed_i = divmod(rem, n_seeds)
    return eps_list[eps_i], nfe_list[nfe_i], seed_list[seed_i]


def run_cell(eps, nfe, seed, metric_cap, output_dir):
    np.random.seed(seed)

    def evaluate(unit_dvs):
        dvs = DV_LO + unit_dvs * (DV_HI - DV_LO)
        return analytic_objectives(dvs, ANTI_IDEAL).tolist(), []

    opt = run_optimization(
        algorithm="borg_mm",
        evaluate=evaluate,
        n_dvs=K,
        n_objs=K + 1,
        n_constrs=0,
        epsilons=[eps] * (K + 1),
        nfe=nfe,
        seed=seed,
        output_dir=output_dir,
        # n_islands auto-picked by borg_runner._auto_islands
    )

    if opt.pareto_dvs.shape[0] == 0:
        return {
            "epsilon": eps, "nfe": nfe, "seed": seed,
            "n_solutions": 0, "wall_seconds": opt.elapsed_s,
            "algorithm": opt.algorithm,
            "L2_star_discrepancy": float("nan"),
            "nn_cv": float("nan"), "nn_mean": float("nan"),
            "hyperplane_max_dev": float("nan"),
        }

    physical_dvs = DV_LO + opt.pareto_dvs * (DV_HI - DV_LO)
    objs = opt.pareto_objs

    lb = np.full(K, DV_LO)
    ub = np.full(K, DV_HI)
    if metric_cap > 0 and len(physical_dvs) > metric_cap:
        rng = np.random.default_rng(seed)
        sel = rng.choice(len(physical_dvs), metric_cap, replace=False)
        metric_points = physical_dvs[sel]
    else:
        metric_points = physical_dvs
    metrics = coverage_metrics(metric_points, lb, ub)
    hp_dev = float(np.max(np.abs(objs.sum(axis=1) - ANTI_IDEAL.sum())))
    return {
        "epsilon": eps, "nfe": nfe, "seed": seed,
        "n_solutions": int(len(physical_dvs)),
        "wall_seconds": opt.elapsed_s,
        "algorithm": opt.algorithm,
        "L2_star_discrepancy": float(metrics["L2_star_discrepancy"]),
        "nn_cv": float(metrics["nn_cv"]),
        "nn_mean": float(metrics["nn_mean"]),
        "hyperplane_max_dev": hp_dev,
    }


def aggregate_cells(cells):
    groups: Dict = {}
    for c in cells:
        groups.setdefault((c["epsilon"], c["nfe"]), []).append(c)
    rows = []
    for (eps, nfe), recs in sorted(groups.items()):
        m = lambda f: float(np.nanmean([r[f] for r in recs]))
        s = lambda f: float(np.nanstd([r[f] for r in recs]))
        rows.append({
            "epsilon": eps, "nfe": nfe, "n_seeds": len(recs),
            "n_solutions_mean": m("n_solutions"),
            "n_solutions_std": s("n_solutions"),
            "L2_star_mean": m("L2_star_discrepancy"),
            "L2_star_std": s("L2_star_discrepancy"),
            "nn_cv_mean": m("nn_cv"),
            "nn_cv_std": s("nn_cv"),
            "wall_seconds_mean": m("wall_seconds"),
        })
    return rows


def recommend(aggregated, min_n=400):
    largest_nfe = max(r["nfe"] for r in aggregated)
    candidates = [r for r in aggregated
                  if r["nfe"] == largest_nfe and r["n_solutions_mean"] >= min_n]
    if not candidates:
        candidates = [r for r in aggregated if r["nfe"] == largest_nfe]
    best = min(candidates, key=lambda r: r["nn_cv_mean"])
    return {
        "recommended_epsilon": best["epsilon"],
        "recommended_nfe": best["nfe"],
        "n_solutions_mean": best["n_solutions_mean"],
        "nn_cv_mean": best["nn_cv_mean"],
        "L2_star_mean": best["L2_star_mean"],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--mode", choices=["cell", "aggregate"], required=True)
    p.add_argument("--task-id", type=int, help="Cell index for --mode cell")
    p.add_argument("--epsilons", type=float, nargs="+", default=DEFAULT_EPSILONS)
    p.add_argument("--nfes", type=int, nargs="+", default=DEFAULT_NFES)
    p.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    p.add_argument("--metric-cap", type=int, default=DEFAULT_METRIC_CAP)
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    cells_dir = out_dir / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "algorithm": "borg_mm",
        "epsilons": args.epsilons, "nfes": args.nfes, "seeds": args.seeds,
        "metric_cap": args.metric_cap,
        "n_cells": len(args.epsilons) * len(args.nfes) * len(args.seeds),
    }, indent=2))

    if args.mode == "cell":
        if args.task_id is None:
            sys.exit("--mode cell requires --task-id")
        eps, nfe, seed = _unpack_cell(args.task_id, args.epsilons, args.nfes, args.seeds)
        print(f"[01/eps_nfe_sweep] cell {args.task_id}: eps={eps} nfe={nfe} seed={seed}")
        cell_dir = cells_dir / f"cell_{args.task_id:04d}"
        cell_dir.mkdir(parents=True, exist_ok=True)
        rec = run_cell(eps, nfe, seed, args.metric_cap, cell_dir)
        (cells_dir / f"cell_{args.task_id:04d}.json").write_text(json.dumps(rec, indent=2))
        print(f"  n={rec['n_solutions']} L2*={rec['L2_star_discrepancy']:.5f} "
              f"NN_CV={rec['nn_cv']:.4f} t={rec['wall_seconds']:.1f}s")
        return

    # aggregate
    cells = [json.loads(p.read_text()) for p in sorted(cells_dir.glob("cell_*.json"))]
    if not cells:
        sys.exit(f"No cells found in {cells_dir}")
    agg = aggregate_cells(cells)
    rec = recommend(agg)
    (out_dir / "aggregate.json").write_text(json.dumps({
        "aggregated": agg, "recommendation": rec,
    }, indent=2))
    print("=== Aggregate ===")
    for r in agg:
        print(f"  eps={r['epsilon']:.3f} nfe={r['nfe']:>6d} "
              f"n={r['n_solutions_mean']:7.1f} "
              f"L2*={r['L2_star_mean']:.5f} NN_CV={r['nn_cv_mean']:.4f}")
    print("=== Recommendation ===")
    for k, v in rec.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
