"""Core MOEA-FIND experiment runner.

Split out of the former ``src.experiment_utils`` god module: this file
owns only the optimization orchestration (``run_experiment``). Data
prep lives in :mod:`src.experiment.data`, anti-ideal/epsilon helpers in
:mod:`src.experiment.anti_ideal`, and figure assembly in
:mod:`src.experiment.plots`.
"""

from pathlib import Path
from typing import Optional

import numpy as np

from src.metrics.objectives import (
    compute_ssi_drought_characteristics,
    drought_objectives,
    make_ssi_calculator,
    flows_to_series,
)
from src.discovery.analysis import coverage_metrics
from src.optimization.borg_runner import run_optimization, OptimizationResult
from src.optimization.constraints import (
    ConstraintConfig,
    ConstraintResult,
    compute_all_constraints,
)
from src.optimization.constraints_dv import DVUniformityConfig, compute_dv_constraint
from src.experiment.anti_ideal import build_epsilons


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
    algorithm: str = "borg_mm",
    constraint_cfg: Optional[ConstraintConfig] = None,
    dv_constraint_cfg: Optional[DVUniformityConfig] = None,
    output_dir: Optional[Path] = None,
    chars_fn=None,
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
        algorithm: Backend name (``"borg_mm"`` for production, ``"borg_serial"``
            for single-process runs). Overridden by ``MOEA_FIND_ALGORITHM`` env var.
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
    from src.metrics.drought_metrics import metric_names, resolve_metric_set

    np.random.seed(seed)

    if constraint_cfg is not None and dv_constraint_cfg is not None:
        raise ValueError(
            "constraint_cfg (hydrologic) and dv_constraint_cfg (DV-space) "
            "are mutually exclusive. Pass exactly one."
        )

    metric_set = resolve_metric_set(objective_keys)
    obj_names = metric_names(metric_set)

    # First-event family: hard "≥1 critical SSI-3 event" constraint replaces
    # the legacy 0.0 fallback so empty-event candidates can't game the L1
    # device. Detected by the first_event_ name prefix.
    is_first_event_family = any(
        m.name.startswith("first_event_") for m in metric_set
    )

    n_dvs = generator.n_dvs
    n_objs = len(metric_set) + 1  # +1 for Manhattan norm
    n_constrs = 0
    eval_count = [0]
    infeasible_count = [0]

    # Pre-fit SSI calculator on historical data (skipped when chars_fn is supplied)
    prefitted_ssi = None
    if chars_fn is None:
        prefitted_ssi = make_ssi_calculator(timescale=timescale)
        hist_series = flows_to_series(monthly_1d, start_date="1950-10-01")
        prefitted_ssi.fit(hist_series)

    # Determine number of hard constraints
    if constraint_cfg is not None:
        enabled = constraint_cfg.enabled
        n_constrs = sum(1 for v in enabled.values() if v)
    elif dv_constraint_cfg is not None and dv_constraint_cfg.enabled:
        n_constrs = 1
    if is_first_event_family:
        n_constrs += 1

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

        # Drought characteristics -> objectives (SSI path or short-block path)
        if chars_fn is not None:
            chars = chars_fn(synthetic_1d, synthetic_2d)
        else:
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

        if is_first_event_family:
            no_event = int(chars.get("first_event_present", 0)) == 0
            hard_violations.append(1.0 if no_event else 0.0)
            if no_event:
                infeasible_count[0] += 1

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

        if chars_fn is not None:
            chars = chars_fn(syn_1d, syn_2d)
        else:
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
