"""Script 03 — Epsilon x NFE sensitivity sweep (manuscript §5, Figure 3).

Sweeps the 3D analytic problem across a grid of epsilon values, NFE budgets,
and random seeds. Intended to be run as a SLURM job array (one task per
(eps, nfe, seed) cell) on HPC. The driver supports three modes:

    * Full sweep (serial on one node): --mode sweep
    * Single-cell for SLURM array:     --mode cell --task-id K
    * Aggregation across cell outputs: --mode aggregate

Cell layout:
    index(eps, nfe, seed) = eps_idx * (n_nfe * n_seeds)
                          + nfe_idx * n_seeds
                          + seed_idx

Per-cell outputs live in outputs/exp03_eps_nfe_sweep/cells/cell_<id>.json.
Aggregate output is outputs/exp03_eps_nfe_sweep/aggregate.json.

Run (local full sweep):
    python scripts/03_eps_nfe_sweep.py --mode sweep --metric-cap 1500

Run (SLURM array cell — see 03_eps_nfe_sweep.slurm):
    python scripts/03_eps_nfe_sweep.py --mode cell --task-id $SLURM_ARRAY_TASK_ID

Run (post-array aggregation):
    python scripts/03_eps_nfe_sweep.py --mode aggregate
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from platypus import EpsNSGAII, Problem, Real  # noqa: E402

from src.analysis import coverage_metrics  # noqa: E402
from src.objectives import analytic_objectives  # noqa: E402
from src.plotting.analytic import fig3_eps_nfe_heatmap  # noqa: E402


K = 3
DV_RANGE = (-3.0, 3.0)
ANTI_IDEAL = np.array([3.0, 3.0, 3.0])
OUTPUT_SLUG = "exp03_eps_nfe_sweep"
DEFAULT_EPSILONS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50]
DEFAULT_NFES = [5_000, 20_000, 50_000]
DEFAULT_SEEDS = [42, 17, 91, 5, 13]
DEFAULT_METRIC_CAP = 1500  # O(n^2) discrepancy cap; pass 0 on HPC to disable


def _cell_index(eps_i: int, nfe_i: int, seed_i: int, n_nfe: int, n_seeds: int) -> int:
    return eps_i * n_nfe * n_seeds + nfe_i * n_seeds + seed_i


def _unpack_cell(task_id: int, eps_list, nfe_list, seed_list) -> Tuple[float, int, int]:
    n_nfe = len(nfe_list); n_seeds = len(seed_list)
    eps_i, rem = divmod(task_id, n_nfe * n_seeds)
    nfe_i, seed_i = divmod(rem, n_seeds)
    return eps_list[eps_i], nfe_list[nfe_i], seed_list[seed_i]


def _build_problem() -> Problem:
    def evaluate(variables):
        dvs = np.array([float(v) for v in variables])
        return analytic_objectives(dvs, ANTI_IDEAL).tolist()

    problem = Problem(K, K + 1)
    for i in range(K):
        problem.types[i] = Real(DV_RANGE[0], DV_RANGE[1])
    problem.function = evaluate
    return problem


def run_cell(eps: float, nfe: int, seed: int, metric_cap: int) -> Dict:
    np.random.seed(seed)
    algo = EpsNSGAII(_build_problem(), epsilons=[eps] * (K + 1))
    t0 = time.perf_counter()
    algo.run(nfe)
    wall = time.perf_counter() - t0

    dvs = np.array([s.variables[:] for s in algo.result])
    objs = np.array([s.objectives[:] for s in algo.result])

    lb = np.full(K, DV_RANGE[0])
    ub = np.full(K, DV_RANGE[1])
    if metric_cap > 0 and len(dvs) > metric_cap:
        rng = np.random.default_rng(seed)
        sel = rng.choice(len(dvs), metric_cap, replace=False)
        metric_points = dvs[sel]
    else:
        metric_points = dvs
    metrics = coverage_metrics(metric_points, lb, ub)

    hp_dev = float(np.max(np.abs(objs.sum(axis=1) - ANTI_IDEAL.sum())))

    return {
        "epsilon": eps,
        "nfe": nfe,
        "seed": seed,
        "n_solutions": int(len(dvs)),
        "wall_seconds": wall,
        "L2_star_discrepancy": float(metrics["L2_star_discrepancy"]),
        "nn_cv": float(metrics["nn_cv"]),
        "nn_mean": float(metrics["nn_mean"]),
        "hyperplane_max_dev": hp_dev,
    }


def aggregate_cells(cells: List[Dict]) -> List[Dict]:
    groups: Dict = {}
    for c in cells:
        groups.setdefault((c["epsilon"], c["nfe"]), []).append(c)

    rows = []
    for (eps, nfe), recs in sorted(groups.items()):
        def m(f): return float(np.mean([r[f] for r in recs]))
        def s(f): return float(np.std([r[f] for r in recs]))
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


def recommend(aggregated: List[Dict], min_n: int = 400) -> Dict:
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
    p.add_argument("--mode", choices=["sweep", "cell", "aggregate"], default="sweep")
    p.add_argument("--epsilons", type=float, nargs="+", default=DEFAULT_EPSILONS)
    p.add_argument("--nfes", type=int, nargs="+", default=DEFAULT_NFES)
    p.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    p.add_argument("--task-id", type=int, help="Cell index (for --mode cell)")
    p.add_argument("--metric-cap", type=int, default=DEFAULT_METRIC_CAP,
                   help="L2* subsample cap; 0 = no cap (use on HPC)")
    p.add_argument("--output-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--plot", action="store_true",
                   help="After aggregation, write manuscript Figure 3 to figures/.")
    p.add_argument("--figure-path", type=Path,
                   default=PROJECT_ROOT / "figures" / "fig03_eps_nfe_sweep.pdf")
    args = p.parse_args()

    cells_dir = args.output_dir / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)

    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "03_eps_nfe_sweep.py",
        "manuscript_section": "§5 Analytic Validation (Fig 3)",
        "epsilons": args.epsilons, "nfes": args.nfes, "seeds": args.seeds,
        "metric_cap": args.metric_cap,
    }, indent=2))

    if args.mode == "cell":
        if args.task_id is None:
            sys.exit("--mode cell requires --task-id")
        eps, nfe, seed = _unpack_cell(args.task_id, args.epsilons, args.nfes, args.seeds)
        print(f"[03] cell {args.task_id}: eps={eps} nfe={nfe} seed={seed}")
        rec = run_cell(eps, nfe, seed, args.metric_cap)
        (cells_dir / f"cell_{args.task_id:04d}.json").write_text(
            json.dumps(rec, indent=2))
        print(f"  n={rec['n_solutions']} "
              f"L2*={rec['L2_star_discrepancy']:.5f} "
              f"NN_CV={rec['nn_cv']:.4f} "
              f"t={rec['wall_seconds']:.1f}s")
        return

    if args.mode == "aggregate":
        cells = [json.loads(p.read_text()) for p in sorted(cells_dir.glob("cell_*.json"))]
        if not cells:
            sys.exit(f"No cells found in {cells_dir}")
        agg = aggregate_cells(cells)
        rec = recommend(agg)
        (args.output_dir / "aggregate.json").write_text(json.dumps({
            "aggregated": agg, "recommendation": rec,
        }, indent=2))
        print("=== Aggregate ===")
        for r in agg:
            print(f"  eps={r['epsilon']:.3f} nfe={r['nfe']:>6d} "
                  f"n={r['n_solutions_mean']:7.1f} "
                  f"L2*={r['L2_star_mean']:.5f} "
                  f"NN_CV={r['nn_cv_mean']:.4f}")
        print("=== Recommendation ===")
        for k, v in rec.items():
            print(f"  {k}: {v}")
        if args.plot:
            import matplotlib.pyplot as plt
            fig = fig3_eps_nfe_heatmap(agg)
            args.figure_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(args.figure_path, dpi=300)
            plt.close(fig)
            print(f"[03] figure: {args.figure_path}")
        return

    # Full sweep (local, single-process)
    total = len(args.epsilons) * len(args.nfes) * len(args.seeds)
    print(f"[03] sweep: {total} cells (eps x nfe x seed)")
    cells: List[Dict] = []
    for ei, eps in enumerate(args.epsilons):
        for ni, nfe in enumerate(args.nfes):
            for si, seed in enumerate(args.seeds):
                cid = _cell_index(ei, ni, si, len(args.nfes), len(args.seeds))
                rec = run_cell(eps, nfe, seed, args.metric_cap)
                (cells_dir / f"cell_{cid:04d}.json").write_text(json.dumps(rec, indent=2))
                cells.append(rec)
                print(f"  cell {cid:3d}: eps={eps:.3f} nfe={nfe} seed={seed} "
                      f"n={rec['n_solutions']} L2*={rec['L2_star_discrepancy']:.5f} "
                      f"NN_CV={rec['nn_cv']:.4f} t={rec['wall_seconds']:.1f}s")
    agg = aggregate_cells(cells)
    rec = recommend(agg)
    (args.output_dir / "aggregate.json").write_text(json.dumps({
        "aggregated": agg, "recommendation": rec,
    }, indent=2))
    print("=== Recommendation ===")
    for k, v in rec.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
