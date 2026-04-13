"""Script 02 — Analytic 3D proof-of-concept (manuscript §5, Figure 2).

Runs the 3-objective Manhattan-norm analytic problem (k=3, 4 objectives
including the Manhattan norm) and verifies that the Pareto front tiles the
2-simplex in objective space. Produces the hyperplane-check and coverage-
comparison numbers cited in DD-10.

Outputs under outputs/exp02_analytic_3d/:
    - results.json
    - pareto.npz
    - config.json

Run:
    python scripts/02_analytic_3d.py --nfe 50000 --seed 42 --plot
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

from platypus import EpsNSGAII, Problem, Real  # noqa: E402

from src.analysis import (  # noqa: E402
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)
from src.objectives import analytic_objectives  # noqa: E402
from src.plotting.analytic import fig2_3d_projections  # noqa: E402


K = 3
DV_RANGE = (-3.0, 3.0)
ANTI_IDEAL = np.array([3.0, 3.0, 3.0])
DEFAULT_EPSILON = 0.15
OUTPUT_SLUG = "exp02_analytic_3d"


def run(nfe: int, seed: int, epsilon: float) -> dict:
    np.random.seed(seed)

    def evaluate(variables):
        dvs = np.array([float(v) for v in variables])
        return analytic_objectives(dvs, ANTI_IDEAL).tolist()

    problem = Problem(K, K + 1)
    for i in range(K):
        problem.types[i] = Real(DV_RANGE[0], DV_RANGE[1])
    problem.function = evaluate

    algo = EpsNSGAII(problem, epsilons=[epsilon] * (K + 1))
    t0 = time.perf_counter()
    algo.run(nfe)
    wall = time.perf_counter() - t0

    dvs = np.array([s.variables[:] for s in algo.result])
    objs = np.array([s.objectives[:] for s in algo.result])

    lb = np.full(K, DV_RANGE[0])
    ub = np.full(K, DV_RANGE[1])
    n = len(dvs)
    lhs = generate_lhs_samples(n, K, lb, ub, seed=seed)
    sobol = generate_sobol_samples(n, K, lb, ub, seed=seed)[:n]
    rng = np.random.default_rng(seed)
    rand = lb + rng.random((n, K)) * (ub - lb)

    obj_sums = objs.sum(axis=1)
    expected_sum = float(ANTI_IDEAL.sum())

    return {
        "n_solutions": int(n),
        "wall_seconds": wall,
        "hyperplane_max_dev": float(np.max(np.abs(obj_sums - expected_sum))),
        "coverage": {
            "pareto": coverage_metrics(dvs, lb, ub),
            "lhs": coverage_metrics(lhs, lb, ub),
            "sobol": coverage_metrics(sobol, lb, ub),
            "random": coverage_metrics(rand, lb, ub),
        },
        "pareto_dvs": dvs,
        "pareto_objs": objs,
    }


def plot(pareto_dvs: np.ndarray, out_path: Path) -> None:
    """Produce the manuscript Figure 2 panel-b via src.plotting.analytic."""
    import matplotlib.pyplot as plt

    fig = fig2_3d_projections(pareto_dvs)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--nfe", type=int, default=50_000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    p.add_argument("--plot", action="store_true")
    p.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / OUTPUT_SLUG)
    p.add_argument("--figure-path", type=Path,
                   default=PROJECT_ROOT / "figures" / "fig02_analytic_3d.pdf")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps({
        "script": "02_analytic_3d.py",
        "manuscript_section": "§5 Analytic Validation (Fig 2, DD-10)",
        "nfe": args.nfe, "seed": args.seed, "epsilon": args.epsilon,
    }, indent=2))

    print(f"[02] 3D analytic: nfe={args.nfe} seed={args.seed} eps={args.epsilon}")
    r = run(args.nfe, args.seed, args.epsilon)

    summary = {k: v for k, v in r.items() if k not in ("pareto_dvs", "pareto_objs")}
    (args.output_dir / "results.json").write_text(json.dumps(summary, indent=2, default=float))
    np.savez(args.output_dir / "pareto.npz", dvs=r["pareto_dvs"], objs=r["pareto_objs"])

    cov = r["coverage"]
    print(f"  n={r['n_solutions']} "
          f"hp_dev={r['hyperplane_max_dev']:.1e}")
    print(f"  L2* pareto={cov['pareto']['L2_star_discrepancy']:.5f} "
          f"lhs={cov['lhs']['L2_star_discrepancy']:.5f} "
          f"sobol={cov['sobol']['L2_star_discrepancy']:.5f}")
    print(f"  NN_CV pareto={cov['pareto']['nn_cv']:.4f} "
          f"lhs={cov['lhs']['nn_cv']:.4f} "
          f"sobol={cov['sobol']['nn_cv']:.4f}")

    if args.plot:
        plot(r["pareto_dvs"], args.figure_path)
        print(f"  figure: {args.figure_path}")


if __name__ == "__main__":
    main()
