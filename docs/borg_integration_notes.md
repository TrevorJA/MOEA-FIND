# Borg MOEA Integration Notes

Findings from debugging Borg MM integration on Hopper HPC (April 2026).

## borg.py Version

Use `~/training/MMBorgBootcamp/Step1-DTLZ2-Optimization/borg.py` — standard borg.py with `seed=` kwarg in `Borg.__init__`.

Do NOT use `~/Research/MMBorgMOEA/borg.py` (passNFE_ALH_PyCheckpoint branch). That branch has a 6-arg C callback `(v, o, c, NFE, island, rank)` but the pre-compiled `libborgmm.so` expects the standard 3-arg `(v, o, c)`. Mismatched callback causes segfaults during `solveMPI`.

## Callback Signature

```python
def func(*vars):
    return (list_of_objectives, list_of_constraints)  # must be a tuple
```

- DVs unpacked as positional args (standard borg.py)
- Return a **tuple**, not a flat list — flat list is interpreted as objectives-only

## Library Loading

`borg.py` loads `./libborg.so` relative to CWD. Workaround in `borg_runner.py`: `os.chdir(lib/borg/)` before import and `startMPI()`, then restore.

## MPI Launch

Per Water Programming bootcamp (`~/training/MMBorgBootcamp/Step1-DTLZ2-Optimization/parallel_multiseeds.sh`):

```bash
mpirun --oversubscribe -np $n_processors python script.py
```

- `--oversubscribe` is required
- Do NOT use `srun --mpi=pmi2` — Borg calls `MPI_Init` internally which conflicts
- `Configuration.startMPI()` must be called BEFORE `Borg(...)` constructor

## SLURM Directives

```bash
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
```

## Bugs Fixed in borg_runner.py

1. Wrapper returned flat list instead of tuple — crashes with constraints
2. `_parse_borg_result` used `sol.getObjectives()` to extract constraints — should use `sol.getConstraints()`
3. `Borg()` created before `startMPI()` — C problem object belonged to wrong library
4. `len(result)` in `_parse_borg_result` — Borg's `Result` class has `size()` not `__len__`.
   `len()` raised `TypeError`, silently dropping ALL Pareto solutions. Fix: `result.size()`.

## Multi-Node MPI

Hopper's Ethernet (eno1np0, 128.84.3.x) is firewalled between compute nodes.
MPI's OOB runtime layer defaults to Ethernet and fails with TCP connection errors.
Fix: route through InfiniBand (ib0, 192.168.12.x):
```bash
mpirun --oversubscribe -np N \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python script.py
```

Do NOT `module purge` — default modules include networking components
needed for multi-node MPI. Just `module load python/3.11.5`.

## Month Alignment (calendar-year vs water-year)

SynHydro's KirschGenerator outputs flows in **calendar-year order** (Jan=col 0,
Dec=col 11). The rest of MOEA-FIND uses **water-year order** (Oct=col 0,
Sep=col 11) because `prepare_data()` aligns historical flows to Oct-Sep.

Fix in `kirsch_wrapper.py`: `np.roll(synthetic, 3, axis=1)` after generation
converts calendar→water year. Without this, the seasonal cycle plot is shifted
~3 months.

## Other Fixes

- `_common.sh`: module names (`python/3.11.5`, `openmpi4/4.0.5`), `gnu9/9.3.0` dependency, account/partition, `CLUSTER_VENV` ordering
- All SLURM scripts: `$(dirname "$0")` → `${SLURM_SUBMIT_DIR}/workflows/` (spool copy issue)
- `constraints.py` + `diag_constraint_calibration.py`: SSI alignment bug (`series.loc[ssi.index]`)
- `pyproject.toml`: invalid build backend

## MM Borg topology and NFE semantics (resolved 2026-04-28)

Two separate findings forced a wrapper redesign while porting the
analytic stage to MM Borg:

1. **Always run with ``n_islands >= 2``.** With ``n_islands == 1``
   (master-slave mode under the multi-master library), Borg
   intermittently SIGSEGVs on a worker rank during MPI startup once
   the worker count exceeds ~6. Reproducible across multiple seeds and
   compute nodes. The fix is in
   :func:`src.borg_runner._auto_islands` -- it returns
   ``max(2, n_ranks // 16)`` for any allocation of >= 5 ranks.
   Drivers should pass ``n_islands=None`` to opt into the heuristic;
   the production stage 04 slurm overrides with ``--n-islands 4`` per
   the Hadka & Reed 2015 production-scale ratio.

2. **``solveMPI(maxEvaluations=...)`` is per-island, not total.** Per
   the Water Programming Python-wrapper tutorial, ``maxEvaluations``
   is the budget delivered to each island master, so the total NFE
   across the run is ``islands * maxEvaluations``. Our wrapper now
   takes the caller's ``nfe`` as the total budget and divides by
   ``n_islands`` before passing to ``solveMPI`` so the call site has
   the obvious meaning. The legacy code was silently giving the
   stage-04 production runs ``4x`` the requested NFE.

## Concurrent-mpirun stress (open issue)

Empirically, when many SLURM array tasks each launch their own
``mpirun -np 8`` on the same node concurrently (e.g., 36 cells x 8
ranks = 288 OpenMPI processes hitting one node simultaneously), Borg
MM segfaults with elevated frequency: ~13/36 cells in a single
unthrottled batch on Hopper, with failures clustered on whichever
nodes ended up host to the most simultaneous mpiruns. Mitigation:
throttle the SLURM array via ``#SBATCH --array=0-N%M`` so concurrent
mpirun count stays bounded (e.g., ``%8`` keeps it to 8 cells x 8
ranks = 64 ranks per dispatch). The throttling is documented in the
stage 01 slurm files.

## Variant Slug System

Output directories are keyed by a deterministic slug encoding all experiment
parameters: `{mode}_T{n_years}_nfe{nfe}_s{seed}_{constrained|unconstrained}`.
Built by `make_variant_slug()` in `src/experiment_utils.py`. Figures saved
per-variant in `{variant_dir}/figures/` — no collisions between configs.
