"""Analytic 2D proof-of-concept (manuscript Fig 1).

Runs the 2-objective Manhattan-norm analytic problem under MM Borg MOEA
(launched via MPI) and writes numerical artifacts only. Figures are
produced separately by ``workflows/01_analytic_validation/plots/analytic_2d.py``.

Outputs under ``outputs/01_analytic_validation/analytic_2d/<slug>/``:
    - config.json  (invocation record)
    - results.json (summary, coverage metrics)
    - pareto.npz   (DVs, objectives)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from src.analysis import (  # noqa: E402
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)
from src.borg_runner import run_optimization  # noqa: E402
from src.objectives import analytic_objectives  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "analytic_2d"
K = 2
DV_LO, DV_HI = -3.0, 3.0
ANTI_IDEAL = np.array([3.0, 3.0])


def slug(nfe: int, seed: int, epsilon: float) -> str:
    return f"k2_nfe{nfe}_eps{epsilon:.3f}_s{seed}"


def run(nfe: int, seed: int, epsilon: float, output_dir: Path) -> dict:
    np.random.seed(seed)

    def evaluate(unit_dvs: np.ndarray):
        # Borg DVs live in [0,1]; scale to the physical [-3,3] domain.
        dvs = DV_LO + unit_dvs * (DV_HI - DV_LO)
        objs = analytic_objectives(dvs, ANTI_IDEAL).tolist()
        return objs, []

    opt = run_optimization(
        algorithm="borg_mm",
        evaluate=evaluate,
        n_dvs=K,
        n_objs=K + 1,
        n_constrs=0,
        epsilons=[epsilon] * (K + 1),
        nfe=nfe,
        seed=seed,
        output_dir=output_dir,
        # n_islands auto-picked by borg_runner._auto_islands
    )

    if opt.pareto_dvs.shape[0] == 0:
        return {
            "n_solutions": 0,
            "wall_seconds": opt.elapsed_s,
            "algorithm": opt.algorithm,
            "pareto_dvs": np.empty((0, K)),
            "pareto_objs": np.empty((0, K + 1)),
        }

    physical_dvs = DV_LO + opt.pareto_dvs * (DV_HI - DV_LO)
    objs = opt.pareto_objs

    lb = np.full(K, DV_LO)
    ub = np.full(K, DV_HI)
    n = len(physical_dvs)
    lhs = generate_lhs_samples(n, K, lb, ub, seed=seed)
    sobol = generate_sobol_samples(n, K, lb, ub, seed=seed)[:n]

    obj_sums = objs.sum(axis=1)
    expected_sum = float(ANTI_IDEAL.sum())

    return {
        "n_solutions": int(n),
        "wall_seconds": opt.elapsed_s,
        "algorithm": opt.algorithm,
        "hyperplane_max_dev": float(np.max(np.abs(obj_sums - expected_sum))),
        "coverage": {
            "pareto": coverage_metrics(physical_dvs, lb, ub),
            "lhs": coverage_metrics(lhs, lb, ub),
            "sobol": coverage_metrics(sobol, lb, ub),
        },
        "pareto_dvs": physical_dvs,
        "pareto_objs": objs,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--epsilon", type=float, required=True)
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER, slug(args.nfe, args.seed, args.epsilon))
    (out_dir / "config.json").write_text(json.dumps({
        "stage": STAGE, "driver": DRIVER,
        "algorithm": "borg_mm",
        "nfe": args.nfe, "seed": args.seed, "epsilon": args.epsilon,
        "K": K, "dv_range": [DV_LO, DV_HI], "anti_ideal": ANTI_IDEAL.tolist(),
    }, indent=2))

    print(f"[01/analytic_2d] nfe={args.nfe} seed={args.seed} eps={args.epsilon}")
    r = run(args.nfe, args.seed, args.epsilon, out_dir)

    if r["n_solutions"] == 0:
        print("  WARNING: no Pareto solutions; nothing to write.")
        return

    summary = {k: v for k, v in r.items() if k not in ("pareto_dvs", "pareto_objs")}
    (out_dir / "results.json").write_text(json.dumps(summary, indent=2, default=float))
    np.savez(out_dir / "pareto.npz", dvs=r["pareto_dvs"], objs=r["pareto_objs"])

    print(f"  n={r['n_solutions']} "
          f"hp_dev={r['hyperplane_max_dev']:.1e} "
          f"L2*={r['coverage']['pareto']['L2_star_discrepancy']:.5f} "
          f"t={r['wall_seconds']:.1f}s")
    print(f"  outputs: {out_dir}")


if __name__ == "__main__":
    main()
