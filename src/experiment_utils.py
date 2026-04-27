"""Shared utilities for MOEA-FIND Kirsch experiment scripts.

Provides common data loading, SSI computation, MOEA execution, and plotting
functions used by both proof-of-concept and kirsch_ensemble experiment scripts.

The ``run_experiment`` function is the primary entry point. It:
  1. Builds an evaluate callback that maps DVs -> (objectives, constraints).
  2. Delegates to :func:`src.borg_runner.run_optimization` for algorithm dispatch
     (MM Borg, serial Borg, or EpsNSGAII dev fallback).
  3. Post-processes the result into a JSON-serialisable dict.
"""

import time
from pathlib import Path
from typing import Dict, Optional, Tuple

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
from src.borg_runner import run_optimization, OptimizationResult
from src.constraints import ConstraintConfig, ConstraintResult, compute_all_constraints
from src.constraints_dv import DVUniformityConfig, compute_dv_constraint


# ---------------------------------------------------------------------------
# Data loading helpers (unchanged)
# ---------------------------------------------------------------------------

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
    objective_keys,
    headroom: float = 1.5,
    feasible_maxes: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """Build the DD-11 anti-ideal point ``D*`` from historical drought chars.

    Thin shim that resolves ``objective_keys`` to a tuple of
    :class:`src.drought_metrics.DroughtMetric` and delegates to
    :func:`src.drought_metrics.compute_anti_ideal`. Each metric's
    :class:`AntiIdealRule` decides placement: ``HEADROOM_TIMES_MAX`` for
    unbounded-above non-cyclic metrics, ``CYCLIC_HEADROOM`` for cyclic
    calendar metrics (``12 × headroom``), and ``CONSTANT`` for metrics
    with a natural upper bound (e.g. fractions in ``[0, 1]``).

    DD-11 requires ``D_j(x) <= D*_j`` for every feasible ``x``; the
    rules above guarantee this provided ``headroom > 1``.

    Args:
        hist_chars: Historical drought characteristics dict.
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`DroughtMetric` instances.
        headroom: Safety factor applied to the historical maximum.
        feasible_maxes: Optional ``{metric_name: observed_max}`` override
            for ``HEADROOM_TIMES_MAX`` metrics. Cyclic and constant
            metrics ignore this argument.

    Returns:
        Anti-ideal array ``D*`` of length ``len(objective_keys)``.
    """
    from src.drought_metrics import compute_anti_ideal, resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    return compute_anti_ideal(
        metric_set,
        hist_chars,
        headroom=headroom,
        feasible_maxes=feasible_maxes,
    )


def extract_pareto_maxes(
    results_json_path: Path,
    objective_keys,
) -> Dict[str, float]:
    """Read a prior ``results.json`` and return ``{name: max_D_j}`` per objective.

    Used to drive a Pareto-based anti-ideal placement on subsequent runs
    via :func:`compute_ssi_anti_ideal`'s ``feasible_maxes`` argument.

    Args:
        results_json_path: Path to a prior ``results.json``.
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.drought_metrics.DroughtMetric` instances.

    Raises:
        ValueError: ``objective_keys`` don't match the reference file's
            objective ordering. Mismatched objective sets would silently
            mis-scale ``D*``.
    """
    import json
    from src.drought_metrics import metric_names, resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    names = metric_names(metric_set)

    payload = json.loads(Path(results_json_path).read_text())
    ref_keys = tuple(payload.get("objective_keys", ()))
    if names != ref_keys:
        raise ValueError(
            f"objective_keys mismatch: requested {names} but "
            f"{results_json_path} was run with {ref_keys}. Rebuild the "
            f"reference run with matching objectives, or drop the --anti-ideal-"
            f"reference flag to fall back to historical max."
        )
    dm = np.asarray(payload.get("drought_metrics", []), dtype=float)
    if dm.size == 0:
        raise ValueError(
            f"{results_json_path} has no drought_metrics (empty Pareto); "
            f"cannot derive feasible maxes."
        )
    return {n: float(dm[:, j].max()) for j, n in enumerate(names)}


# ---------------------------------------------------------------------------
# Epsilon defaults
# ---------------------------------------------------------------------------


def build_epsilons(
    objective_keys,
    epsilon_map: Optional[Dict[str, float]] = None,
    manhattan_eps: Optional[float] = None,
) -> list:
    """Return epsilon list for objectives plus the Manhattan-norm auxiliary.

    Args:
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.drought_metrics.DroughtMetric` instances. When
            metric instances are passed, the per-axis epsilon is read
            from each metric's ``epsilon`` field. When string names are
            passed, ``epsilon_map`` is consulted (falls back to ``0.5``).
        epsilon_map: Optional override mapping metric name → epsilon.
            Used only when ``objective_keys`` is a tuple of strings.
        manhattan_eps: Epsilon for the ``f_{K+1}`` auxiliary objective.
            Defaults to
            :data:`src.experiment_config.DEFAULT_EXPERIMENT.manhattan_eps`.
    """
    from src.drought_metrics import DroughtMetric, resolve_metric_set

    keys = tuple(objective_keys)
    if len(keys) > 0 and isinstance(keys[0], DroughtMetric):
        eps = [m.epsilon for m in keys]
    elif epsilon_map is not None:
        eps = [epsilon_map.get(k, 0.5) for k in keys]
    else:
        # Resolve string names to metrics so each carries its own epsilon.
        metric_set = resolve_metric_set(keys)
        eps = [m.epsilon for m in metric_set]

    if manhattan_eps is None:
        from src.experiment_config import DEFAULT_EXPERIMENT
        manhattan_eps = DEFAULT_EXPERIMENT.manhattan_eps
    eps.append(float(manhattan_eps))
    return eps


# ---------------------------------------------------------------------------
# Variant identification
# ---------------------------------------------------------------------------

def make_variant_slug(
    mode: str,
    n_years: int,
    nfe: int,
    seed: int,
    constrained: bool,
    *,
    extra: Optional[dict] = None,
) -> str:
    """Build a deterministic, filesystem-safe variant identifier.

    The slug encodes all parameters that distinguish experiment runs so that
    each unique configuration gets its own output directory.

    Args:
        mode: DV injection mode (``"index"`` or ``"residual"``).
        n_years: Synthetic trace length in years.
        nfe: Maximum function evaluations.
        seed: Random seed.
        constrained: Whether plausibility constraints are active.
        extra: Optional additional knobs, e.g. ``{"ssi": "6"}``.
            Keys are sorted for determinism.

    Returns:
        Slug like ``"residual_T20_nfe500000_s42_constrained"``.
    """
    parts = [
        mode,
        f"T{n_years}",
        f"nfe{nfe}",
        f"s{seed}",
        "constrained" if constrained else "unconstrained",
    ]
    if extra:
        for k, v in sorted(extra.items()):
            parts.append(f"{k}{v}")
    return "_".join(parts)


# ---------------------------------------------------------------------------
# Core experiment runner
# ---------------------------------------------------------------------------

def run_experiment(
    monthly_2d: np.ndarray,
    monthly_1d: np.ndarray,
    generator,
    generator_name: str,
    objective_keys,
    anti_ideal: np.ndarray,
    timescale: int,
    n_years_out: int,
    nfe: int,
    seed: int,
    algorithm: str = "eps_nsga2",
    constraint_cfg: Optional[ConstraintConfig] = None,
    dv_constraint_cfg: Optional[DVUniformityConfig] = None,
    output_dir: Optional[Path] = None,
    **algo_kwargs,
) -> dict:
    """Run a single MOEA experiment with SSI objectives and plausibility constraints.

    This is the main entry point for all Kirsch-based MOEA-FIND experiments.
    It builds an evaluation callback that maps decision variables to
    (objectives, constraints), delegates to :func:`borg_runner.run_optimization`,
    then post-processes into a serialisable result dict.

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
        algorithm: Backend name (``"borg_mm"``, ``"borg_serial"``,
            ``"eps_nsga2"``). Overridden by ``MOEA_FIND_ALGORITHM`` env var.
        constraint_cfg: Optional :class:`ConstraintConfig`. If None,
            hydrologic constraints are disabled.
        dv_constraint_cfg: Optional :class:`DVUniformityConfig` for the
            DV-space uniformity ablation arm. Mutually exclusive with
            ``constraint_cfg`` — setting both raises ValueError. If None,
            DV-space constraint is disabled.
        output_dir: Directory for algorithm runtime/checkpoint files.
            Defaults to a temp directory.
        **algo_kwargs: Forwarded to the backend (e.g. ``n_islands``).

    Returns:
        Results dict with Pareto front, hyperplane check, ranges,
        coverage metrics, and constraint diagnostics.
    """
    from src.drought_metrics import metric_names, resolve_metric_set

    np.random.seed(seed)

    if constraint_cfg is not None and dv_constraint_cfg is not None:
        raise ValueError(
            "constraint_cfg (hydrologic) and dv_constraint_cfg (DV-space) "
            "are mutually exclusive. Pass exactly one."
        )

    metric_set = resolve_metric_set(objective_keys)
    obj_names = metric_names(metric_set)

    n_dvs = generator.n_dvs
    n_objs = len(metric_set) + 1  # +1 for Manhattan norm
    n_constrs = 0
    eval_count = [0]
    infeasible_count = [0]

    # Pre-fit SSI calculator on historical data
    prefitted_ssi = make_ssi_calculator(timescale=timescale)
    hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
    prefitted_ssi.fit(hist_series)

    # Determine number of hard constraints
    if constraint_cfg is not None:
        enabled = constraint_cfg.enabled
        n_constrs = sum(1 for v in enabled.values() if v)
    elif dv_constraint_cfg is not None and dv_constraint_cfg.enabled:
        n_constrs = 1

    epsilons = build_epsilons(metric_set)

    if output_dir is None:
        import tempfile
        output_dir = Path(tempfile.mkdtemp(prefix="moea_find_"))

    def evaluate(dvs: np.ndarray):
        """Callback: DVs -> (objectives, constraints).

        Returns objectives with soft penalty folded into Manhattan norm,
        and hard constraints in Borg's <= 0 convention.
        """
        eval_count[0] += 1
        synthetic_2d = generator.generate(dvs)
        if synthetic_2d.ndim == 1:
            synthetic_2d = synthetic_2d.reshape(n_years_out, 12)
        synthetic_1d = synthetic_2d.flatten()

        # SSI -> drought characteristics -> objectives
        syn_series = flows_to_series(synthetic_1d, start_date="2100-01-01")
        ssi_syn = prefitted_ssi.transform(syn_series)
        chars = compute_ssi_drought_characteristics(
            ssi_syn, monthly_flows=synthetic_1d
        )
        objs = list(drought_objectives(chars, anti_ideal, metric_set))

        # Constraints
        hard_violations = []
        if constraint_cfg is not None:
            cr = compute_all_constraints(
                synthetic_1d, synthetic_2d, constraint_cfg,
                ssi_calc=prefitted_ssi,
            )
            hard_violations = cr.hard_violations
            if not cr.feasible:
                infeasible_count[0] += 1
            # Fold soft penalty into Manhattan norm (last objective)
            objs[-1] += cr.soft_penalty_weighted
        elif dv_constraint_cfg is not None:
            cr_dv = compute_dv_constraint(dvs, dv_constraint_cfg)
            hard_violations = cr_dv.hard_violations
            if not cr_dv.feasible:
                infeasible_count[0] += 1
            objs[-1] += cr_dv.soft_penalty_weighted

        return objs, hard_violations

    # --- Run optimisation ---
    print(f"  Running {algorithm} ({generator_name}): "
          f"{n_dvs} DVs, {n_objs} obj, {n_constrs} constr, {nfe} NFE ...")

    opt_result: OptimizationResult = run_optimization(
        algorithm=algorithm,
        evaluate=evaluate,
        n_dvs=n_dvs,
        n_objs=n_objs,
        n_constrs=n_constrs,
        epsilons=epsilons,
        nfe=nfe,
        seed=seed,
        output_dir=output_dir,
        **algo_kwargs,
    )

    elapsed = opt_result.elapsed_s
    print(f"  Done in {elapsed:.1f}s ({eval_count[0]} evals, "
          f"{elapsed / max(eval_count[0], 1) * 1000:.1f} ms/eval)")
    if constraint_cfg is not None or dv_constraint_cfg is not None:
        print(f"  Infeasible evals: {infeasible_count[0]} / {eval_count[0]} "
              f"({100 * infeasible_count[0] / max(eval_count[0], 1):.1f}%)")

    # --- Post-process ---
    pareto_objs = opt_result.pareto_objs
    pareto_dvs = opt_result.pareto_dvs
    n_pareto = len(pareto_objs)
    print(f"  Pareto solutions: {n_pareto}")

    if n_pareto == 0:
        return {
            "mode": generator_name,
            "algorithm": opt_result.algorithm,
            "ssi_timescale": timescale,
            "n_dvs": n_dvs,
            "n_years_out": n_years_out,
            "nfe": nfe,
            "elapsed_s": elapsed,
            "n_pareto": 0,
            "n_evals_total": eval_count[0],
            "n_infeasible": infeasible_count[0],
            "anti_ideal": anti_ideal.tolist(),
            "epsilons": epsilons,
            "objective_keys": list(obj_names),
            "drought_metrics": [],
            "pareto_chars": [],
            "pareto_dvs": [],
            "pareto_traces_1d": [],
            "pareto_traces_2d": [],
            "error": "No Pareto-optimal solutions found",
        }

    drought_metrics = pareto_objs[:, :len(metric_set)]

    # Hyperplane check
    obj_sums = np.sum(pareto_objs, axis=1)
    expected_sum = np.sum(anti_ideal)
    print(f"  Hyperplane: expected={expected_sum:.2f}, "
          f"mean={np.mean(obj_sums):.4f}, std={np.std(obj_sums):.6f}")

    for j, name in enumerate(obj_names):
        print(f"  {name}: [{drought_metrics[:, j].min():.2f}, "
              f"{drought_metrics[:, j].max():.2f}]")

    lb = np.zeros(len(metric_set))
    ub = anti_ideal.copy()
    dm = coverage_metrics(drought_metrics, lb, ub)

    # Regenerate Pareto traces for diagnostics and plotting
    pareto_chars = []
    pareto_traces_1d = []
    pareto_traces_2d = []
    pareto_constraint_diags = []
    for dvs_row in pareto_dvs:
        syn_2d = generator.generate(dvs_row)
        if syn_2d.ndim == 1:
            syn_2d = syn_2d.reshape(n_years_out, 12)
        syn_1d = syn_2d.flatten()

        pareto_traces_1d.append(syn_1d.tolist())
        pareto_traces_2d.append(syn_2d.tolist())

        syn_s = flows_to_series(syn_1d, start_date="2100-01-01")
        ssi_re = prefitted_ssi.transform(syn_s)
        chars = compute_ssi_drought_characteristics(
            ssi_re, monthly_flows=syn_1d
        )
        pareto_chars.append(chars)

        if constraint_cfg is not None:
            cr = compute_all_constraints(
                syn_1d, syn_2d, constraint_cfg, ssi_calc=prefitted_ssi
            )
            pareto_constraint_diags.append(cr.to_dict())
        elif dv_constraint_cfg is not None:
            cr_dv = compute_dv_constraint(dvs_row, dv_constraint_cfg)
            pareto_constraint_diags.append(cr_dv.to_dict())

    result = {
        "mode": generator_name,
        "algorithm": opt_result.algorithm,
        "ssi_timescale": timescale,
        "n_dvs": n_dvs,
        "n_years_out": n_years_out,
        "nfe": nfe,
        "elapsed_s": elapsed,
        "n_pareto": n_pareto,
        "n_evals_total": eval_count[0],
        "n_infeasible": infeasible_count[0],
        "anti_ideal": anti_ideal.tolist(),
        "epsilons": epsilons,
        "objective_keys": list(obj_names),
        "hyperplane": {
            "expected_sum": float(expected_sum),
            "actual_mean": float(np.mean(obj_sums)),
            "actual_std": float(np.std(obj_sums)),
        },
        "ranges": {
            name: {
                "min": float(drought_metrics[:, j].min()),
                "max": float(drought_metrics[:, j].max()),
            }
            for j, name in enumerate(obj_names)
        },
        "coverage": dm,
        "drought_metrics": drought_metrics.tolist(),
        "pareto_chars": pareto_chars,
        "pareto_dvs": pareto_dvs.tolist(),
        "pareto_traces_1d": pareto_traces_1d,
        "pareto_traces_2d": pareto_traces_2d,
    }
    if constraint_cfg is not None:
        result["constraint_mode"] = "hydrologic"
        result["constraint_config"] = constraint_cfg.to_dict()
        result["constraint_diagnostics"] = pareto_constraint_diags
    elif dv_constraint_cfg is not None:
        result["constraint_mode"] = "dv_uniform"
        result["dv_constraint_config"] = dv_constraint_cfg.to_dict()
        result["constraint_diagnostics"] = pareto_constraint_diags
    else:
        result["constraint_mode"] = "none"

    return result


# ---------------------------------------------------------------------------
# Plotting (unchanged from previous version)
# ---------------------------------------------------------------------------

def plot_comparison(
    results_list: list,
    hist_chars: dict,
    anti_ideal: np.ndarray,
    objective_keys,
    fig_dir: Path,
) -> None:
    """Generate comparison figures for multiple generators.

    Creates scatter plots of Pareto fronts in drought space, overlaying
    historical point and anti-ideal point.

    Args:
        results_list: List of results dicts from run_experiment.
        hist_chars: Historical drought characteristics dict.
        anti_ideal: Anti-ideal point (k-dimensional).
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.drought_metrics.DroughtMetric` instances.
        fig_dir: Output directory for figures.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.drought_metrics import resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    if len(metric_set) < 2:
        raise ValueError("plot_comparison requires at least two metrics")

    m0, m1 = metric_set[0], metric_set[1]

    fig_dir.mkdir(parents=True, exist_ok=True)

    color_map = {
        "Kirsch (index)": "#d62728",
        "Kirsch (residual)": "#9467bd",
    }

    n_modes = len(results_list)
    fig, axes = plt.subplots(1, n_modes, figsize=(6 * n_modes, 5), squeeze=False)
    axes = axes[0]

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
            m0.extract(hist_chars), m1.extract(hist_chars),
            marker="*", s=200, c="black", zorder=5, label="Historical",
        )
        ax.scatter(
            anti_ideal[0], anti_ideal[1],
            marker="x", s=200, c="red", zorder=5, label="Anti-ideal D*",
        )
        ax.set_xlabel(f"{m0.label} ({m0.units})")
        ax.set_ylabel(f"{m1.label} ({m1.units})")
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
          f"{m0.label[:16]:>16} {m1.label[:16]:>16} {'L2*':>8}")
    for r in results_list:
        rng = r["ranges"]
        cov = r["coverage"]
        l2 = cov.get("L2_star_discrepancy", cov.get("L2_star", 0))
        r0 = rng.get(m0.name, {"min": 0.0, "max": 0.0})
        r1 = rng.get(m1.name, {"min": 0.0, "max": 0.0})
        print(f"  {r['mode']:<25} {r['n_pareto']:>5} "
              f"{r0['min']:>7.2f}-{r0['max']:<7.2f} "
              f"{r1['min']:>7.2f}-{r1['max']:<7.2f} "
              f"{l2:>8.4f}")
