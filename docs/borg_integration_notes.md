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

## Variant Slug System

Output directories are keyed by a deterministic slug encoding all experiment
parameters: `{mode}_T{n_years}_nfe{nfe}_s{seed}_{constrained|unconstrained}`.
Built by `make_variant_slug()` in `src/experiment_utils.py`. Figures saved
per-variant in `{variant_dir}/figures/` — no collisions between configs.
