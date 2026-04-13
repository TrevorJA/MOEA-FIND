"""Shared utilities for MOEA-FIND Kirsch experiment scripts.

Provides common data loading, SSI computation, MOEA execution, and plotting
functions used by both proof-of-concept and kirsch_ensemble experiment scripts.
"""

import time
from pathlib import Path
from typing import Tuple

import numpy as np

from src.data import load_usgs_daily, daily_to_monthly
from src.objectives import (
    compute_ssi,
    compute_ssi_drought_characteristics,
    drought_objectives,
    make_ssi_calculator,
    flows_to_series,
)
from src.analysis import coverage_metrics


def prepare_data(cache_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load historical USGS data and return (monthly_2d, monthly_1d).

    Args:
        cache_dir: Directory for caching downloaded USGS data.

    Returns:
        Tuple of (monthly_2d, monthly_1d) where:
        - monthly_2d: shape (n_years, 12)
        - monthly_1d: flattened historical flows
    """
    daily = load_usgs_daily(cache_dir=cache_dir)
    monthly = daily_to_monthly(daily)

    # Align to water year (October-September)
    first_oct = monthly.index[monthly.index.month == 10][0]
    last_sep = monthly.index[monthly.index.month == 9][-1]
    monthly = monthly[first_oct:last_sep]

    n_years = len(monthly) // 12
    monthly_values = monthly.values[:n_years * 12]
    monthly_2d = monthly_values.reshape(n_years, 12)

    print(f"Historical: {n_years} water years, mean={monthly_2d.mean():.1f} cfs")
    return monthly_2d, monthly_values


def compute_historical_ssi_chars(
    monthly_1d,
    timescale: int,
    **kwargs,
) -> Tuple:
    """Compute SSI and drought characteristics for historical data.

    Args:
        monthly_1d: 1D array or Series of monthly flows.
        timescale: SSI accumulation period (1, 3, 6, 12).

    Returns:
        Tuple of (ssi_series, ssi_calculator, drought_characteristics_dict).
    """
    ssi, calc = compute_ssi(monthly_1d, timescale=timescale)
    chars = compute_ssi_drought_characteristics(ssi)
    return ssi, calc, chars


def compute_ssi_anti_ideal(
    hist_chars: dict,
    objective_keys: tuple,
    headroom: float = 1.5,
) -> np.ndarray:
    """Compute anti-ideal from historical SSI drought characteristics.

    Sets anti-ideal at headroom * max observed value for each metric.

    Args:
        hist_chars: Historical drought characteristics dict from
            compute_ssi_drought_characteristics.
        objective_keys: Metric keys for objectives.
        headroom: Multiplier for anti-ideal placement (default 1.5).

    Returns:
        Anti-ideal array (k-dimensional).
    """
    max_key_map = {
        "mean_duration": "max_duration",
        "mean_magnitude": "max_magnitude",
        "mean_severity": "worst_severity",
        "mean_avg_severity": "worst_severity",
        "max_duration": "max_duration",
        "max_magnitude": "max_magnitude",
        "worst_severity": "worst_severity",
        "frequency": "frequency",
    }

    anti_ideal = []
    for key in objective_keys:
        max_key = max_key_map.get(key, key)
        val = hist_chars.get(max_key, 10.0)
        if val == 0:
            val = 10.0  # fallback for zero-event case
        anti_ideal.append(val * headroom)

    return np.array(anti_ideal)


def run_experiment(
    monthly_2d: np.ndarray,
    monthly_1d: np.ndarray,
    generator,
    generator_name: str,
    objective_keys: tuple,
    anti_ideal: np.ndarray,
    timescale: int,
    n_years_out: int,
    nfe: int,
    seed: int,
) -> dict:
    """Run a single MOEA experiment with SSI objectives.

    Args:
        monthly_2d: Historical flows (n_years, 12).
        monthly_1d: Historical flows 1D.
        generator: Initialized KirschBorgWrapper instance.
        generator_name: Human-readable name for logging.
        objective_keys: Which SSI drought metrics to optimize.
        anti_ideal: Anti-ideal point (k-dimensional).
        timescale: SSI accumulation period.
        n_years_out: Synthetic trace length in years.
        nfe: Number of function evaluations.
        seed: Random seed.

    Returns:
        Results dict with Pareto front, hyperplane check, ranges, coverage metrics.
    """
    from platypus import EpsNSGAII, Problem, Real

    np.random.seed(seed)

    # KirschBorgWrapper has n_dvs as a property
    n_dvs = generator.n_dvs
    n_objs = len(objective_keys) + 1  # +1 for Manhattan norm
    eval_count = [0]

    # Pre-fit SSI calculator on historical data
    prefitted_ssi = make_ssi_calculator(timescale=timescale)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    prefitted_ssi.fit(hist_series)

    def evaluate(variables):
        """Objective function for Borg/EpsNSGAII.

        Takes decision variables, generates synthetic flows, computes SSI and
        drought metrics, returns objectives (metrics + Manhattan norm).
        """
        eval_count[0] += 1
        dvs = np.array([float(v) for v in variables])

        # Generate synthetic flows via KirschBorgWrapper
        synthetic_2d = generator.generate(dvs)

        # Ensure 2D shape (n_years_out, 12)
        if synthetic_2d.ndim == 1:
            synthetic_2d = synthetic_2d.reshape(n_years_out, 12)

        # Transform to SSI using pre-fitted calculator
        syn_series = flows_to_series(
            synthetic_2d.flatten(), start_date="2100-01-01",
        )
        ssi_syn = prefitted_ssi.transform(syn_series)
        chars = compute_ssi_drought_characteristics(ssi_syn)

        # Return objectives: metrics + Manhattan norm
        return list(drought_objectives(chars, anti_ideal, objective_keys))

    # Set up Borg problem
    problem = Problem(n_dvs, n_objs)
    for i in range(n_dvs):
        problem.types[i] = Real(0.0, 1.0)
    problem.function = evaluate

    # Epsilon values (grid spacing for Borg's epsilon-dominance)
    eps_map = {
        "mean_duration": 0.5,
        "mean_magnitude": 0.5,
        "mean_severity": 0.1,
        "mean_avg_severity": 0.1,
        "max_duration": 1.0,
        "max_magnitude": 1.0,
        "worst_severity": 0.1,
        "frequency": 0.5,
    }
    epsilons = [eps_map.get(k, 0.5) for k in objective_keys]
    eps_manhattan = sum(epsilons)
    epsilons.append(eps_manhattan)

    # Run EpsNSGAII
    algorithm = EpsNSGAII(problem, epsilons=epsilons)

    print(f"  Running EpsNSGAII ({generator_name}): "
          f"{n_dvs} DVs, {n_objs} obj, {nfe} NFE...")
    t0 = time.time()
    algorithm.run(nfe)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s "
          f"({eval_count[0]} evals, {elapsed/max(eval_count[0],1)*1000:.1f} ms/eval)")

    # Extract results
    pareto_objs = np.array([list(s.objectives) for s in algorithm.result])
    pareto_dvs = np.array([
        [float(v) for v in s.variables] for s in algorithm.result
    ])
    print(f"  Pareto solutions: {len(pareto_objs)}")

    if len(pareto_objs) == 0:
        return {"n_pareto": 0, "mode": generator_name, "error": "No solutions"}

    drought_metrics = pareto_objs[:, :len(objective_keys)]

    # Hyperplane check (Manhattan norm trick should force solutions to hyperplane)
    obj_sums = np.sum(pareto_objs, axis=1)
    expected_sum = np.sum(anti_ideal)
    print(f"  Hyperplane: expected={expected_sum:.2f}, "
          f"mean={np.mean(obj_sums):.4f}, std={np.std(obj_sums):.6f}")

    # Objective ranges
    for j, key in enumerate(objective_keys):
        print(f"  {key}: [{drought_metrics[:, j].min():.2f}, "
              f"{drought_metrics[:, j].max():.2f}]")

    # Coverage metrics in normalized space
    lb = np.zeros(len(objective_keys))
    ub = anti_ideal.copy()
    dm = coverage_metrics(drought_metrics, lb, ub)

    # Regenerate Pareto traces for full characteristics (for diagnostics)
    pareto_chars = []
    for dvs in pareto_dvs:
        syn_2d = generator.generate(dvs)
        if syn_2d.ndim == 1:
            syn_2d = syn_2d.reshape(n_years_out, 12)
        syn_s = flows_to_series(syn_2d.flatten(), start_date="2100-01-01")
        ssi_re = prefitted_ssi.transform(syn_s)
        chars = compute_ssi_drought_characteristics(ssi_re)
        pareto_chars.append(chars)

    return {
        "mode": generator_name,
        "ssi_timescale": timescale,
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
        "ranges": {
            key: {
                "min": float(drought_metrics[:, j].min()),
                "max": float(drought_metrics[:, j].max()),
            }
            for j, key in enumerate(objective_keys)
        },
        "coverage": dm,
        "drought_metrics": drought_metrics.tolist(),
        "pareto_chars": pareto_chars,
    }


def plot_comparison(
    results_list: list,
    hist_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys: tuple,
    fig_dir: Path,
) -> None:
    """Generate comparison figures for multiple generators.

    Creates scatter plots of Pareto fronts in drought space, overlaying
    historical point and anti-ideal point.

    Args:
        results_list: List of results dicts from run_experiment.
        hist_chars: Historical drought characteristics dict.
        anti_ideal: Anti-ideal point (k-dimensional).
        objective_keys: Objective names for axis labels.
        fig_dir: Output directory for figures.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)

    color_map = {
        "Kirsch (index)": "#d62728",
        "Kirsch (residual)": "#9467bd",
    }

    n_modes = len(results_list)
    fig, axes = plt.subplots(1, n_modes, figsize=(6 * n_modes, 5), squeeze=False)
    axes = axes[0]

    k0, k1 = objective_keys[0], objective_keys[1]

    for ax, r in zip(axes, results_list):
        metrics = np.array(r["drought_metrics"])
        mode = r["mode"]
        c = color_map.get(mode, "gray")

        ax.scatter(
            metrics[:, 0], metrics[:, 1],
            s=15, alpha=0.7, c=c,
            label=f"{mode} (n={r['n_pareto']})",
        )
        ax.scatter(
            hist_chars[k0], hist_chars[k1],
            marker="*", s=200, c="black", zorder=5, label="Historical",
        )
        ax.scatter(
            anti_ideal[0], anti_ideal[1],
            marker="x", s=200, c="red", zorder=5, label="Anti-ideal D*",
        )
        ax.set_xlabel(f"{k0} (months)")
        ax.set_ylabel(f"{k1} (avg severity)")
        ax.set_title(f"{mode} (n={r['n_pareto']})")
        ax.legend(fontsize=8)

    ssi_acc = results_list[0].get("ssi_timescale", "?")
    fig.suptitle(f"Kirsch: SSI-{ssi_acc} Drought Space (index vs residual)", fontsize=13)
    fig.tight_layout()

    fname = f"kirsch_poc_ssi{ssi_acc}_coverage.png"
    fig.savefig(fig_dir / fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname}")

    # Summary table
    print(f"\n  === SSI-{ssi_acc} Kirsch PoC Summary ===")
    print(f"  {'Generator':<25} {'N':>5} "
          f"{'Duration range':>16} {'Severity range':>16} {'L2*':>8}")
    for r in results_list:
        rng = r["ranges"]
        cov = r["coverage"]
        l2 = cov.get("L2_star_discrepancy", cov.get("L2_star", 0))
        print(f"  {r['mode']:<25} {r['n_pareto']:>5} "
              f"{rng[k0]['min']:.1f}-{rng[k0]['max']:.1f} mo      "
              f"{rng[k1]['min']:.3f}-{rng[k1]['max']:.3f}     "
              f"{l2:>8.4f}")
