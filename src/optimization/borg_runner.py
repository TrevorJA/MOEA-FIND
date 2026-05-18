"""Algorithm runners for MOEA-FIND: MM Borg (production) and serial Borg.

Provides a single entry point :func:`run_optimization` that dispatches to the
requested Borg backend. Both backends accept the same evaluation callback
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

MM Borg is the production algorithm used for every optimization run in
the manuscript -- analytic, calibration, and Cannonsville. Serial Borg
is the single-process fallback for cases where launching MPI is not
warranted (e.g. a one-off interactive smoke test). No other MOEA backend
is supported.
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

def _auto_islands(n_ranks: int) -> int:
    """Pick a safe ``n_islands`` for ``n_ranks`` total MPI processes.

    With ``n_islands == 1`` (master-slave mode), Borg MM gives one of
    the ranks the controller role and the rest become workers under a
    single master. Empirically this topology SIGSEGVs on a worker rank
    during MPI startup whenever the worker count is large enough to
    make the master/controller queue a bottleneck (reproducible at
    n_ranks >= ~7 across multiple seeds). The fix is to **always run
    multi-master** (``n_islands >= 2``) on any allocation that supports
    it, so the controller has at least two island-masters to talk to
    and per-master worker queues stay short.

    Heuristic:
        - n_ranks < 5  : 1 island (degenerate; only valid for tiny smoke tests)
        - 5 <= n_ranks : ``max(2, n_ranks // 16)`` -- 2 islands at the
          floor, scaling up at ~16 ranks per island (matches the
          production stage 04 ratio of 4 islands at 120 ranks).

    The MPI rank requirement for ``n_islands`` islands and ``W``
    workers per island is ``n_islands * (W + 1) + 1``; this rule keeps
    each island master with a manageable W = (n_ranks - 1) // n_islands
    - 1 worker count.
    """
    if n_ranks < 5:
        return 1
    return max(2, n_ranks // 16)


def run_borg_mm(
    evaluate: EvalFn,
    n_dvs: int,
    n_objs: int,
    n_constrs: int,
    epsilons: List[float],
    nfe: int,
    seed: int,
    output_dir: Path,
    n_islands: Optional[int] = None,
    checkpoint_freq: int = 10000,
    old_checkpoint: Optional[str] = None,
) -> OptimizationResult:
    """Run multi-master Borg MOEA via MPI.

    Must be launched under ``mpirun`` / ``srun``. ``nfe`` is the
    **total** function-evaluation budget across all islands; this
    wrapper divides by ``n_islands`` internally because Borg's
    ``solveMPI(maxEvaluations=...)`` is per-island.

    Args:
        evaluate: Evaluation callback. Takes a 1-D DV array; returns
            (objectives, constraints) where constraints follow the
            Borg ``<= 0`` feasible convention.
        n_dvs: Number of decision variables.
        n_objs: Number of objectives.
        n_constrs: Number of constraints.
        epsilons: Epsilon values per objective.
        nfe: Total function-evaluation budget (split equally across islands).
        seed: Random seed.
        output_dir: Directory for runtime files and checkpoints.
        n_islands: Number of islands. ``None`` (default) lets
            :func:`_auto_islands` pick a topology robust to the
            single-island Borg MM bug.
        checkpoint_freq: NFE interval between checkpoint writes.
        old_checkpoint: Path to a checkpoint file to resume from.
    """
    # Resolve rank count from the MPI environment to inform the islands
    # heuristic. OpenMPI sets OMPI_COMM_WORLD_SIZE; SLURM sets SLURM_NTASKS.
    n_ranks = int(
        os.environ.get("OMPI_COMM_WORLD_SIZE")
        or os.environ.get("SLURM_NTASKS")
        or "1"
    )
    if n_islands is None:
        n_islands = _auto_islands(n_ranks)

    # borg.py loads ./libborg.so relative to CWD, so cd there for import.
    _borg_dir = str(Path(__file__).resolve().parents[2] / "lib" / "borg")
    _saved_cwd = os.getcwd()
    os.chdir(_borg_dir)
    from borg import Borg, Configuration  # type: ignore[import-not-found]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Standard borg.py unpacks DVs as positional args; return (objs, constrs) tuple.
    def _mm_wrapper(*vars_args):
        dvs = np.array(vars_args, dtype=float)
        objs, constrs = evaluate(dvs)
        return (list(objs), list(constrs))

    # startMPI() loads ./libborgmm.so from CWD, so stay in _borg_dir.
    Configuration.startMPI()
    os.chdir(_saved_cwd)

    borg = Borg(
        numberOfVariables=n_dvs,
        numberOfObjectives=n_objs,
        numberOfConstraints=n_constrs,
        function=_mm_wrapper,
        epsilons=epsilons,
        bounds=[[0.0, 1.0]] * n_dvs,
        seed=seed,
    )

    # Borg's solveMPI(maxEvaluations=...) is the per-island budget;
    # divide by n_islands so the caller's ``nfe`` is the total.
    per_island_nfe = max(1, nfe // max(1, n_islands))
    print(f"[borg_runner] mm topology: n_ranks={n_ranks} n_islands={n_islands} "
          f"per_island_nfe={per_island_nfe} (total≈{per_island_nfe * n_islands})",
          flush=True)

    t0 = time.time()
    result = borg.solveMPI(
        islands=n_islands,
        maxEvaluations=per_island_nfe,
        frequency=checkpoint_freq,
    )
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
    **_ignored,
) -> OptimizationResult:
    """Run serial (single-process) Borg MOEA.

    No MPI required. Useful for local debugging with the Borg C library
    but without a cluster allocation.
    """
    _borg_dir = str(Path(__file__).resolve().parents[2] / "lib" / "borg")
    _saved_cwd = os.getcwd()
    os.chdir(_borg_dir)
    from borg import Borg, Configuration  # type: ignore[import-not-found]
    os.chdir(_saved_cwd)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Standard borg.py unpacks DVs as positional args; return (objs, constrs) tuple.
    def _serial_wrapper(*vars_args):
        dvs = np.array(vars_args, dtype=float)
        objs, constrs = evaluate(dvs)
        return (list(objs), list(constrs))

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
    if result is None or result.size() == 0:
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
        objs_list.append(sol.getObjectives())
        if n_constrs > 0:
            constrs_list.append(sol.getConstraints())
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
# Dispatcher
# ---------------------------------------------------------------------------

ALGORITHM_CHOICES = ("borg_mm", "borg_serial")


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
    raise RuntimeError("unreachable")
