"""MOEA-FIND Proof-of-Concept: Analytic test problems in 2D and 3D.

Experiment 1.1: 2D analytic (X1, X2, Manhattan) — validates Manhattan norm trick
Experiment 1.2: 3D analytic (X1, X2, X3, Manhattan) — tests uniformity on 2-simplex

Uses platypus EpsNSGAII as a local stand-in for Borg MOEA.
Both use epsilon-dominance archiving, the key mechanism for uniform tiling.

Usage:
    python experiments/proof_of_concept/run_analytic_poc.py [--nfe 100000] [--seed 42]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from platypus import (
    EpsNSGAII,
    Problem,
    Real,
)

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.objectives import analytic_objectives, normalize_to_unit_cube
from src.analysis import (
    coverage_metrics,
    generate_lhs_samples,
    generate_sobol_samples,
)


def run_experiment(
    k: int,
    dv_range: Tuple[float, float],
    anti_ideal: np.ndarray,
    epsilons: List[float],
    nfe: int,
    seed: int,
) -> Dict:
    """Run analytic POC for k-dimensional problem.

    Args:
        k: Number of drought-like objectives (DVs = objectives 1..k).
        dv_range: (lower, upper) bounds for each DV.
        anti_ideal: Anti-ideal point (k-dimensional).
        epsilons: Epsilon values for each of k+1 objectives.
        nfe: Number of function evaluations.
        seed: Random seed.

    Returns:
        Dict with Pareto front points, metrics, and comparison data.
    """
    np.random.seed(seed)

    # Define platypus problem
    # platypus calls function(variables) and expects list of objectives returned
    def evaluate(variables):
        dvs = np.array([float(v) for v in variables])
        objs = analytic_objectives(dvs, anti_ideal)
        return objs.tolist()

    problem = Problem(k, k + 1)
    for i in range(k):
        problem.types[i] = Real(dv_range[0], dv_range[1])
    problem.function = evaluate

    # Run EpsNSGAII (epsilon-dominance, like Borg)
    algorithm = EpsNSGAII(problem, epsilons=epsilons)
    algorithm.run(nfe)

    # Extract Pareto front
    # DV space = first k objectives (since J_i = X_i)
    pareto_dvs = np.array([s.variables[:] for s in algorithm.result])
    pareto_objs = np.array([s.objectives[:] for s in algorithm.result])

    # Bounds for coverage metrics
    lb = np.full(k, dv_range[0])
    ub = np.full(k, dv_range[1])

    # Compute coverage metrics for Pareto front projected to DV space
    pareto_metrics = coverage_metrics(pareto_dvs, lb, ub)

    # Generate comparison samples with same number of points
    n_pareto = len(pareto_dvs)
    lhs_samples = generate_lhs_samples(n_pareto, k, lb, ub, seed=seed)
    lhs_metrics = coverage_metrics(lhs_samples, lb, ub)

    sobol_samples = generate_sobol_samples(n_pareto, k, lb, ub, seed=seed)
    # Sobol may return more points (power of 2); trim to match
    sobol_samples = sobol_samples[:n_pareto]
    sobol_metrics = coverage_metrics(sobol_samples, lb, ub)

    # Random uniform baseline
    rng = np.random.default_rng(seed)
    random_samples = lb + rng.random((n_pareto, k)) * (ub - lb)
    random_metrics = coverage_metrics(random_samples, lb, ub)

    # Verify hyperplane property: sum of all objectives should be constant
    obj_sums = np.sum(pareto_objs, axis=1)
    expected_sum = np.sum(anti_ideal)

    results = {
        "k": k,
        "n_objectives": k + 1,
        "nfe": nfe,
        "n_pareto_solutions": n_pareto,
        "dv_range": list(dv_range),
        "anti_ideal": anti_ideal.tolist(),
        "epsilons": epsilons,
        "hyperplane_check": {
            "expected_sum": float(expected_sum),
            "actual_mean": float(np.mean(obj_sums)),
            "actual_std": float(np.std(obj_sums)),
            "max_deviation": float(np.max(np.abs(obj_sums - expected_sum))),
        },
        "coverage_comparison": {
            "pareto_front": pareto_metrics,
            "lhs": lhs_metrics,
            "sobol": sobol_metrics,
            "random": random_metrics,
        },
        "pareto_dvs": pareto_dvs.tolist(),
        "pareto_objs": pareto_objs.tolist(),
    }
    return results


def run_2d_experiment(nfe: int, seed: int) -> Dict:
    """Experiment 1.1: 2D analytic test."""
    print(f"=== Experiment 1.1: 2D Analytic (k=2, {nfe} NFE) ===")
    k = 2
    anti_ideal = np.array([3.0, 3.0])
    # Epsilon: controls grid spacing on Pareto front
    # Smaller epsilon = more solutions = finer coverage
    eps = 0.06
    epsilons = [eps] * (k + 1)

    results = run_experiment(
        k=k,
        dv_range=(-3.0, 3.0),
        anti_ideal=anti_ideal,
        epsilons=epsilons,
        nfe=nfe,
        seed=seed,
    )

    print(f"  Pareto solutions: {results['n_pareto_solutions']}")
    hp = results["hyperplane_check"]
    print(f"  Hyperplane check: expected sum={hp['expected_sum']:.2f}, "
          f"actual={hp['actual_mean']:.4f} +/- {hp['actual_std']:.4f}")
    cm = results["coverage_comparison"]
    print(f"  L2* discrepancy — Pareto: {cm['pareto_front']['L2_star_discrepancy']:.6f}, "
          f"LHS: {cm['lhs']['L2_star_discrepancy']:.6f}, "
          f"Sobol: {cm['sobol']['L2_star_discrepancy']:.6f}, "
          f"Random: {cm['random']['L2_star_discrepancy']:.6f}")
    print(f"  NN CV — Pareto: {cm['pareto_front']['nn_cv']:.4f}, "
          f"LHS: {cm['lhs']['nn_cv']:.4f}, "
          f"Sobol: {cm['sobol']['nn_cv']:.4f}, "
          f"Random: {cm['random']['nn_cv']:.4f}")
    return results


def run_3d_experiment(nfe: int, seed: int) -> Dict:
    """Experiment 1.2: 3D analytic test."""
    print(f"\n=== Experiment 1.2: 3D Analytic (k=3, {nfe} NFE) ===")
    k = 3
    anti_ideal = np.array([3.0, 3.0, 3.0])
    eps = 0.15  # Coarser epsilon for 3D (more boxes to fill)
    epsilons = [eps] * (k + 1)

    results = run_experiment(
        k=k,
        dv_range=(-3.0, 3.0),
        anti_ideal=anti_ideal,
        epsilons=epsilons,
        nfe=nfe,
        seed=seed,
    )

    print(f"  Pareto solutions: {results['n_pareto_solutions']}")
    hp = results["hyperplane_check"]
    print(f"  Hyperplane check: expected sum={hp['expected_sum']:.2f}, "
          f"actual={hp['actual_mean']:.4f} +/- {hp['actual_std']:.4f}")
    cm = results["coverage_comparison"]
    print(f"  L2* discrepancy — Pareto: {cm['pareto_front']['L2_star_discrepancy']:.6f}, "
          f"LHS: {cm['lhs']['L2_star_discrepancy']:.6f}, "
          f"Sobol: {cm['sobol']['L2_star_discrepancy']:.6f}, "
          f"Random: {cm['random']['L2_star_discrepancy']:.6f}")
    print(f"  NN CV — Pareto: {cm['pareto_front']['nn_cv']:.4f}, "
          f"LHS: {cm['lhs']['nn_cv']:.4f}, "
          f"Sobol: {cm['sobol']['nn_cv']:.4f}, "
          f"Random: {cm['random']['nn_cv']:.4f}")
    return results


def run_epsilon_sensitivity(nfe: int, seed: int) -> Dict:
    """Experiment 1.3: Epsilon sensitivity for 2D case."""
    print(f"\n=== Experiment 1.3: Epsilon Sensitivity (2D, {nfe} NFE) ===")
    k = 2
    anti_ideal = np.array([3.0, 3.0])

    epsilon_values = [0.02, 0.04, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5]
    sensitivity_results = []

    for eps in epsilon_values:
        epsilons = [eps] * (k + 1)
        res = run_experiment(
            k=k,
            dv_range=(-3.0, 3.0),
            anti_ideal=anti_ideal,
            epsilons=epsilons,
            nfe=nfe,
            seed=seed,
        )
        cm = res["coverage_comparison"]["pareto_front"]
        entry = {
            "epsilon": eps,
            "n_solutions": res["n_pareto_solutions"],
            "L2_star_discrepancy": cm["L2_star_discrepancy"],
            "nn_cv": cm["nn_cv"],
            "nn_mean": cm["nn_mean"],
        }
        sensitivity_results.append(entry)
        print(f"  eps={eps:.3f}: n={entry['n_solutions']}, "
              f"L2*={entry['L2_star_discrepancy']:.6f}, "
              f"NN_CV={entry['nn_cv']:.4f}")

    return {"epsilon_sensitivity": sensitivity_results}


def main():
    parser = argparse.ArgumentParser(description="MOEA-FIND Analytic POC")
    parser.add_argument("--nfe", type=int, default=100000,
                        help="Number of function evaluations")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--output-dir", type=str,
                        default=str(project_root / "outputs" / "poc"),
                        help="Output directory for results")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run experiments
    results_2d = run_2d_experiment(args.nfe, args.seed)
    results_3d = run_3d_experiment(args.nfe, args.seed)
    results_eps = run_epsilon_sensitivity(args.nfe, args.seed)

    # Save results (without large arrays for the summary)
    summary = {
        "experiment_1_1_2d": {
            k: v for k, v in results_2d.items()
            if k not in ("pareto_dvs", "pareto_objs")
        },
        "experiment_1_2_3d": {
            k: v for k, v in results_3d.items()
            if k not in ("pareto_dvs", "pareto_objs")
        },
        "experiment_1_3_sensitivity": results_eps,
    }
    summary_path = output_dir / "poc_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to {summary_path}")

    # Save full Pareto fronts as numpy arrays for plotting
    np.savez(
        output_dir / "pareto_2d.npz",
        dvs=np.array(results_2d["pareto_dvs"]),
        objs=np.array(results_2d["pareto_objs"]),
    )
    np.savez(
        output_dir / "pareto_3d.npz",
        dvs=np.array(results_3d["pareto_dvs"]),
        objs=np.array(results_3d["pareto_objs"]),
    )
    print("Pareto front arrays saved.")

    return summary


if __name__ == "__main__":
    main()
