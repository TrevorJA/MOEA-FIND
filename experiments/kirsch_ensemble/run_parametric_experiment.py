"""MOEA-FIND Experiment 3.1: Parametric CDF + vine copula generator.

Tests the parametric generator (kappa4 marginals + D-vine copula) for
drought characteristic space exploration. Uses threshold-based drought metrics.

Uses USGS gauge 01423000 (West Branch Delaware at Walton, Cannonsville inflow).

Usage:
    python experiments/kirsch_ensemble/run_parametric_experiment.py [--nfe 5000]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from platypus import EpsNSGAII, Problem, Real

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.data import load_usgs_daily, daily_to_monthly
from src.parametric import ParametricGenerator
from src.objectives import (
    compute_drought_characteristics,
    drought_objectives,
    manhattan_norm,
)
from src.constraints import compute_all_constraints
from src.analysis import coverage_metrics


def prepare_historical_data(cache_dir: Path):
    """Load and prepare historical monthly flow data.

    Returns:
        Tuple of (monthly_2d, threshold, hist_chars).
    """
    print("Loading USGS data (gauge 01423000, Cannonsville)...")
    daily = load_usgs_daily(cache_dir=cache_dir)
    monthly = daily_to_monthly(daily)
    print(f"  Daily: {len(daily)} records")
    print(f"  Monthly: {len(monthly)} months")

    # Reshape to water years (Oct-Sep)
    first_oct = monthly.index[monthly.index.month == 10][0]
    last_sep = monthly.index[monthly.index.month == 9][-1]
    monthly = monthly[first_oct:last_sep]
    n_months = len(monthly)
    n_years = n_months // 12
    monthly_values = monthly.values[:n_years * 12]
    monthly_2d = monthly_values.reshape(n_years, 12)
    print(f"  Complete water years: {n_years}")

    threshold = float(np.percentile(monthly_values, 20))
    print(f"  Drought threshold (P20): {threshold:.1f} cfs")

    hist_chars = compute_drought_characteristics(monthly_values, threshold, n_years)
    print(f"  Historical: {hist_chars['n_events']} events, "
          f"freq={hist_chars['frequency']:.1f}/decade, "
          f"dur={hist_chars['mean_duration']:.1f} mo, "
          f"int={hist_chars['mean_intensity']:.1f} cfs")

    return monthly_2d, threshold, hist_chars


def compute_anti_ideal(monthly_2d, threshold, objective_keys, headroom=1.5):
    """Compute anti-ideal from historical max event values."""
    from src.objectives import compute_drought_events
    events = compute_drought_events(monthly_2d.flatten(), threshold)

    if len(events) == 0:
        defaults = {
            "mean_duration": 12.0,
            "mean_intensity": 100.0,
            "mean_avg_severity": 1.0,
        }
        return np.array([defaults.get(k, 10.0) for k in objective_keys])

    anti_ideal = []
    for key in objective_keys:
        if key == "mean_duration":
            val = max(e["duration"] for e in events)
        elif key == "mean_intensity":
            val = max(e["intensity"] for e in events)
        elif key == "mean_avg_severity":
            val = 1.0  # SSI-based severity; fallback for threshold-based experiments
        elif key == "frequency":
            n_years = monthly_2d.shape[0]
            val = len(events) / n_years * 10 * headroom
            anti_ideal.append(val)
            continue
        else:
            val = 10.0
        anti_ideal.append(val * headroom)

    return np.array(anti_ideal)


def run_parametric_experiment(
    monthly_2d: np.ndarray,
    threshold: float,
    hist_chars: dict,
    n_years_out: int,
    nfe: int,
    seed: int,
) -> dict:
    """Run parametric generator through MOEA optimization.

    Args:
        monthly_2d: Historical monthly flows (n_years, 12).
        threshold: Drought threshold.
        hist_chars: Historical drought characteristics.
        n_years_out: Synthetic trace length in years.
        nfe: Number of function evaluations.
        seed: Random seed.

    Returns:
        Results dict with Pareto front and coverage metrics.
    """
    np.random.seed(seed)

    generator = ParametricGenerator(
        monthly_2d, marginal_method="kappa4", max_lag=2,
    )
    n_dvs = generator.n_dvs(n_years_out)

    objective_keys = ("mean_duration", "mean_intensity")
    anti_ideal = compute_anti_ideal(monthly_2d, threshold, objective_keys)
    print(f"  Anti-ideal: dur={anti_ideal[0]:.1f} mo, int={anti_ideal[1]:.1f} cfs")

    n_objs = len(objective_keys) + 1  # +1 Manhattan norm
    eval_count = [0]

    def evaluate(variables):
        eval_count[0] += 1
        cdf_probs = np.array([float(v) for v in variables])

        synthetic_2d = generator.generate(cdf_probs, n_years_out)
        synthetic_1d = synthetic_2d.flatten()

        syn_chars = compute_drought_characteristics(
            synthetic_1d, threshold, n_years_out,
        )
        return list(drought_objectives(syn_chars, anti_ideal, objective_keys))

    problem = Problem(n_dvs, n_objs)
    for i in range(n_dvs):
        problem.types[i] = Real(0.0, 1.0)
    problem.function = evaluate

    epsilons = [0.1, 1.0, 1.1]
    algorithm = EpsNSGAII(problem, epsilons=epsilons)

    print(f"  Running EpsNSGAII (parametric): {n_dvs} DVs, {nfe} NFE...")
    t0 = time.time()
    algorithm.run(nfe)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s ({eval_count[0]} evals, "
          f"{elapsed / max(eval_count[0], 1) * 1000:.1f} ms/eval)")

    pareto_objs = np.array([list(s.objectives) for s in algorithm.result])
    pareto_dvs = np.array([[float(v) for v in s.variables] for s in algorithm.result])
    print(f"  Pareto solutions: {len(pareto_objs)}")

    if len(pareto_objs) == 0:
        return {"n_pareto": 0, "mode": "parametric", "error": "No solutions"}

    drought_metrics = pareto_objs[:, :2]

    # Hyperplane check
    obj_sums = np.sum(pareto_objs, axis=1)
    expected_sum = np.sum(anti_ideal)
    print(f"  Hyperplane: expected={expected_sum:.1f}, "
          f"mean={np.mean(obj_sums):.4f}, std={np.std(obj_sums):.4f}")

    lb_drought = np.array([0.0, 0.0])
    ub_drought = anti_ideal.copy()
    dm = coverage_metrics(drought_metrics, lb_drought, ub_drought)

    print(f"  Duration range: [{drought_metrics[:, 0].min():.2f}, "
          f"{drought_metrics[:, 0].max():.2f}] mo")
    print(f"  Intensity range: [{drought_metrics[:, 1].min():.2f}, "
          f"{drought_metrics[:, 1].max():.2f}] cfs")

    # Regenerate traces for diagnostics
    pareto_chars = []
    for dvs in pareto_dvs:
        syn_2d = generator.generate(dvs, n_years_out)
        syn_1d = syn_2d.flatten()
        chars = compute_drought_characteristics(syn_1d, threshold, n_years_out)
        pareto_chars.append(chars)

    return {
        "mode": "parametric",
        "n_dvs": n_dvs,
        "n_years_out": n_years_out,
        "nfe": nfe,
        "elapsed_s": elapsed,
        "n_pareto": len(pareto_objs),
        "anti_ideal": anti_ideal.tolist(),
        "epsilons": epsilons,
        "objective_keys": list(objective_keys),
        "hyperplane": {
            "expected_sum": float(expected_sum),
            "actual_mean": float(np.mean(obj_sums)),
            "actual_std": float(np.std(obj_sums)),
        },
        "drought_range": {
            "duration_min": float(drought_metrics[:, 0].min()),
            "duration_max": float(drought_metrics[:, 0].max()),
            "intensity_min": float(drought_metrics[:, 1].min()),
            "intensity_max": float(drought_metrics[:, 1].max()),
        },
        "coverage": dm,
        "drought_metrics": drought_metrics.tolist(),
        "pareto_chars": pareto_chars,
    }


def main():
    parser = argparse.ArgumentParser(description="Experiment 3.1: Parametric POC")
    parser.add_argument("--nfe", type=int, default=5000,
                        help="Number of function evaluations")
    parser.add_argument("--n-years", type=int, default=30,
                        help="Length of synthetic traces (years)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    args = parser.parse_args()

    output_dir = project_root / "outputs" / "parametric"
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = project_root / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    monthly_2d, threshold, hist_chars = prepare_historical_data(cache_dir)

    # Parametric generator only
    print("\n=== Experiment 3.1: Parametric generator (kappa4 + D-vine) ===")
    results = run_parametric_experiment(
        monthly_2d, threshold, hist_chars,
        args.n_years, args.nfe, args.seed,
    )

    fname = output_dir / "results_parametric.json"
    with open(fname, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {fname}")


if __name__ == "__main__":
    main()
