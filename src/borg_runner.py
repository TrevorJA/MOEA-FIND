"""Algorithm runners for MOEA-FIND: MM Borg, serial Borg, and EpsNSGAII.

Provides a single entry point :func:`run_optimization` that dispatches to the
requested algorithm backend. All backends accept the same evaluation callback
signature and return a standardised :class:`OptimizationResult`.

Borg API follows the ``passNFE_ALH_PyCheckpoint`` branch of the MMBorgMOEA
repository and the Water Programming blog tutorials (Feb/Aug 2025).

Key Borg conventions:
    - Objective function returns ``[*objectives, *constraints]``.
    - Constraint values ``<= 0`` are feasible.
    - Serial wrapper: ``func(*vars)`` returning the combined list.
    - MM wrapper: ``func(vars, NFE)`` with an extra NFE counter argument.
    - MPI ranks needed: ``n_islands * (workers_per_island + 1) + 1``.
    - Checkpoints: ``newCheckpointFileBase`` / ``oldCheckpointFile`` kwargs.

EpsNSGAII (platypus) is retained as a **dev-only** fallback for local
testing without the Borg C library.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class OptimizationResult:
    """Standardised output from any algorithm backend.

    Attributes:
        pareto_dvs: Decision variable matrix (n_solutions, n_dvs).
        pareto_objs: Objective matrix (n_solutions, n_objs).
        pareto_constrs: Constraint matrix (n_solutions, n_constrs). 0 = feasible.
        n_evals: Actual number of function evaluations performed.
        elapsed_s: Wall-clock time in seconds.
        algorithm: Name string for logging / serialisation.
    """
    pareto_dvs: np.ndarray
    pareto_objs: np.ndarray
    pareto_constrs: np.ndarray
    n_evals: int
    elapsed_s: float
    algorithm: str
    extra: Dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Evaluate callback type
# ---------------------------------------------------------------------------
# All backends expect this signature from the caller:
#     evaluate(dvs: np.ndarray) -> Tuple[List[float], List[float]]
#         dvs: 1-D array of decision variables in [0,1].
#         returns: (objectives, hard_constraint_violations)
#             where constraint value <= 0 is feasible (Borg convention).
EvalFn = Callable[[np.ndarray], Tuple[List[float], List[float]]]


# ---------------------------------------------------------------------------
# MM Borg
# ---------------------------------------------------------------------------

def run_borg_mm(
    evaluate: EvalFn,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    epsilons: List[float],
    nfe: int,
    seed: int,
    output_dir: Path,
    n_islands: int = 1,
    checkpoint_freq: int = 10000,
    old_checkpoint: Optional[str] = None,
) -> OptimizationResult:
    """Run multi-master Borg MOEA via MPI.

    Must be launched under ``mpirun`` / ``srun`` with at least
    ``n_islands * (K+1) + 1`` ranks where K is workers per island
    (determined automatically by Borg from available ranks).

    Args:
        evaluate: Evaluation callback. Takes a 1-D DV array; returns
            (objectives, constraints) where constraints follow the
            Borg ``<= 0`` feasible convention.
        n_dvs: Number of decision variables.
        n_objs: Number of objectives.
        n_constrs: Number of constraints.
        epsilons: Epsilon values per objective.
        nfe: Maximum function evaluations.
        seed: Random seed.
        output_dir: Directory for runtime files and checkpoints.
        n_islands: Number of islands (>=1). Use 1 for master-slave.
        checkpoint_freq: NFE interval between checkpoint writes.
        old_checkpoint: Path to a checkpoint file to resume from.
    """
    from borg import Borg, Configuration  # type: ignore[import-not-found]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # MM Borg passes (vars_list, NFE) to the function.
    def _mm_wrapper(vars_list, nfe_count):
        dvs = np.array(vars_list, dtype=float)
        objs, constrs = evaluate(dvs)
        return list(objs) + list(constrs)

    borg = Borg(
        numberOfVariables=n_dvs,
        numberOfObjectives=n_objs,
        numberOfConstraints=n_constrs,
        function=_mm_wrapper,
        epsilons=epsilons,
        bounds=[[0.0, 1.0]] * n_dvs,
        seed=seed,
    )

    Configuration.startMPI()
    t0 = time.time()

    solve_kwargs = dict(
        islands=n_islands,
        maxEvaluations=nfe,
        runtime=str(output_dir / "runtime_%d.txt"),
        frequency=checkpoint_freq,
        newCheckpointFileBase=str(output_dir / "checkpoint"),
    )
    if old_checkpoint is not None:
        solve_kwargs["oldCheckpointFile"] = str(old_checkpoint)

    result = borg.solveMPI(**solve_kwargs)
    elapsed = time.time() - t0
    Configuration.stopMPI()

    return _parse_borg_result(result, n_dvs, n_objs, n_constrs, nfe, elapsed, "borg_mm")


def run_borg_serial(
    evaluate: EvalFn,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    epsilons: List[float],
    nfe: int,
    seed: int,
    output_dir: Path,
    checkpoint_freq: int = 10000,
    old_checkpoint: Optional[str] = None,
) -> OptimizationResult:
    """Run serial (single-process) Borg MOEA.

    No MPI required. Useful for local debugging with the Borg C library
    but without a cluster allocation.
    """
    from borg import Borg  # type: ignore[import-not-found]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Serial Borg passes *vars (unpacked) to the function.
    def _serial_wrapper(*vars_args):
        dvs = np.array(vars_args, dtype=float)
        objs, constrs = evaluate(dvs)
        return list(objs) + list(constrs)

    borg = Borg(
        numberOfVariables=n_dvs,
        numberOfObjectives=n_objs,
        numberOfConstraints=n_constrs,
        function=_serial_wrapper,
        epsilons=epsilons,
        bounds=[[0.0, 1.0]] * n_dvs,
        seed=seed,
    )

    settings = {
        "maxEvaluations": nfe,
        "runtimefile": str(output_dir / "runtime.txt"),
        "runtimeformat": "borg",
        "frequency": checkpoint_freq,
        "newCheckpointFileBase": str(output_dir / "checkpoint"),
    }
    if old_checkpoint is not None:
        settings["oldCheckpointFile"] = str(old_checkpoint)

    t0 = time.time()
    result = borg.solve(settings=settings)
    elapsed = time.time() - t0

    return _parse_borg_result(result, n_dvs, n_objs, n_constrs, nfe, elapsed, "borg_serial")


def _parse_borg_result(
    result,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    nfe: int,
    elapsed: float,
    algo_name: str,
) -> OptimizationResult:
    """Extract Pareto arrays from a Borg result object."""
    if result is None or len(result) == 0:
        return OptimizationResult(
            pareto_dvs=np.empty((0, n_dvs)),
            pareto_objs=np.empty((0, n_objs)),
            pareto_constrs=np.empty((0, n_constrs)),
            n_evals=nfe,
            elapsed_s=elapsed,
            algorithm=algo_name,
        )

    dvs_list = []
    objs_list = []
    constrs_list = []
    for sol in result:
        dvs_list.append(sol.getVariables())
        combined = sol.getObjectives()
        objs_list.append(combined[:n_objs])
        if n_constrs > 0:
            constrs_list.append(combined[n_objs: n_objs + n_constrs])
        else:
            constrs_list.append([])

    return OptimizationResult(
        pareto_dvs=np.array(dvs_list, dtype=float),
        pareto_objs=np.array(objs_list, dtype=float),
        pareto_constrs=np.array(constrs_list, dtype=float) if n_constrs > 0 else np.empty((len(dvs_list), 0)),
        n_evals=nfe,
        elapsed_s=elapsed,
        algorithm=algo_name,
    )


# ---------------------------------------------------------------------------
# EpsNSGAII (dev-only fallback)
# ---------------------------------------------------------------------------

def run_eps_nsga2(
    evaluate: EvalFn,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    epsilons: List[float],
    nfe: int,
    seed: int,
    output_dir: Path,
) -> OptimizationResult:
    """Run platypus EpsNSGAII. DEV ONLY — not for HPC production runs.

    Constraints are handled via Platypus constraint-domination when
    ``n_constrs > 0``.

    Platypus ``Problem.function`` receives the decoded variable *list*
    (not a Solution) and must return objectives; for constraints we must
    override ``Problem.evaluate(solution)`` instead.
    """
    from platypus import EpsNSGAII, Problem, Real  # type: ignore[import-not-found]

    np.random.seed(seed)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    problem = Problem(n_dvs, n_objs, n_constrs)
    for i in range(n_dvs):
        problem.types[i] = Real(0.0, 1.0)

    if n_constrs > 0:
        # Platypus convention: constraint(x) <= 0 is feasible.
        problem.constraints[:] = "<=0"

    if n_constrs > 0:
        # With constraints, override evaluate(Solution) so we can set
        # both objectives and constraints on the Solution object.
        _orig_evaluate = problem.evaluate

        def _eval_with_constraints(solution):
            dvs = np.array([float(v) for v in solution.variables])
            objs, constrs = evaluate(dvs)
            solution.objectives[:] = objs
            solution.constraints[:] = constrs

        problem.evaluate = _eval_with_constraints
    else:
        # Without constraints, the simple function(variables)->objectives
        # pattern suffices.
        def _simple_evaluate(variables):
            dvs = np.array([float(v) for v in variables])
            objs, _ = evaluate(dvs)
            return objs

        problem.function = _simple_evaluate

    algorithm = EpsNSGAII(problem, epsilons=epsilons)

    t0 = time.time()
    algorithm.run(nfe)
    elapsed = time.time() - t0

    solutions = list(algorithm.result)
    if n_constrs > 0:
        feasible = [s for s in solutions if s.feasible]
        if feasible:
            solutions = feasible

    if not solutions:
        return OptimizationResult(
            pareto_dvs=np.empty((0, n_dvs)),
            pareto_objs=np.empty((0, n_objs)),
            pareto_constrs=np.empty((0, n_constrs)),
            n_evals=nfe,
            elapsed_s=elapsed,
            algorithm="eps_nsga2",
        )

    dvs_arr = np.array([[float(v) for v in s.variables] for s in solutions])
    objs_arr = np.array([list(s.objectives) for s in solutions])
    if n_constrs > 0:
        constrs_arr = np.array([list(s.constraints) for s in solutions])
    else:
        constrs_arr = np.empty((len(solutions), 0))

    return OptimizationResult(
        pareto_dvs=dvs_arr,
        pareto_objs=objs_arr,
        pareto_constrs=constrs_arr,
        n_evals=nfe,
        elapsed_s=elapsed,
        algorithm="eps_nsga2",
    )


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

ALGORITHM_CHOICES = ("borg_mm", "borg_serial", "eps_nsga2")


def run_optimization(
    algorithm: str,
    evaluate: EvalFn,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    epsilons: List[float],
    nfe: int,
    seed: int,
    output_dir: Path,
    **kwargs,
) -> OptimizationResult:
    """Unified dispatcher for all supported MOEA backends.

    The ``MOEA_FIND_ALGORITHM`` environment variable overrides ``algorithm``
    if set, so slurm scripts can select the backend without code changes.

    Args:
        algorithm: One of ``ALGORITHM_CHOICES``.
        evaluate: Evaluation callback ``(dvs) -> (objs, constrs)``.
        n_dvs: Number of decision variables.
        n_objs: Number of objectives.
        n_constrs: Number of constraints.
        epsilons: Per-objective epsilon values.
        nfe: Maximum function evaluations.
        seed: Random seed.
        output_dir: Directory for runtime/checkpoint/result files.
        **kwargs: Forwarded to the backend (e.g. ``n_islands``,
            ``checkpoint_freq``, ``old_checkpoint``).
    """
    env_algo = os.environ.get("MOEA_FIND_ALGORITHM", "").strip()
    if env_algo:
        algorithm = env_algo

    if algorithm not in ALGORITHM_CHOICES:
        raise ValueError(
            f"Unknown algorithm {algorithm!r}. Choose from {ALGORITHM_CHOICES}."
        )

    print(f"[borg_runner] algorithm={algorithm}, n_dvs={n_dvs}, n_objs={n_objs}, "
          f"n_constrs={n_constrs}, nfe={nfe}, seed={seed}")

    if algorithm == "borg_mm":
        return run_borg_mm(
            evaluate, n_dvs, n_objs, n_constrs, epsilons, nfe, seed,
            output_dir, **kwargs,
        )
    elif algorithm == "borg_serial":
        return run_borg_serial(
            evaluate, n_dvs, n_objs, n_constrs, epsilons, nfe, seed,
            output_dir, **kwargs,
        )
    elif algorithm == "eps_nsga2":
        return run_eps_nsga2(
            evaluate, n_dvs, n_objs, n_constrs, epsilons, nfe, seed,
            output_dir,
        )
    raise RuntimeError("unreachable")
