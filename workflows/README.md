# MOEA-FIND Workflows

All numerical experiments for the manuscript live in this directory, organized as:

```
workflows/
  _common.sh              # shared SLURM/bash helpers (modules, venv, MPI launcher)
  experiments/             # numbered Python experiment drivers
    01_analytic_2d.py
    ...
  slurm/                   # matching SLURM batch scripts
    01_analytic_2d.slurm
    ...
    slurm_logs/            # SLURM stdout/stderr logs (gitignored)
  diagnostics/             # standalone diagnostic/calibration scripts
    diag_constraint_calibration.py
    ...
```

Each experiment is one numbered Python driver (in `experiments/`) plus one matching SLURM batch script (in `slurm/`). Every entry below is pinned to a specific section and figure of the manuscript draft.

## Manuscript mapping

| ID | Script | Manuscript section | Figure(s) | Parallel model |
|----|--------|--------------------|-----------|----------------|
| 01 | `01_analytic_2d.py` | §5 Analytic Validation | Fig 1, 2 | serial |
| 02 | `02_analytic_3d.py` | §5 Analytic Validation | Fig 2 | serial |
| 03 | `03_eps_nfe_sweep.py` | §5 Analytic Validation | Fig 3 | SLURM job array (one task per (ε, NFE, seed) cell) |
| 04 | `04_kirsch_single_site.py` | §6.1–6.2 Single-Site Kirsch | Fig 5, 6 | MPI master-worker (MM Borg) |
| 05 | `05_kirsch_library_build.py` | §6.3 Library Baseline | Fig 7 | mpi4py work queue |
| 06 | `06_library_subsample_baseline.py` | §6.3 Library Baseline | Fig 7 | serial (vectorized) |
| 07 | `07_event_level_kirsch.py` | §6.4 Event-level Formulation | Fig 5 (inset) | MPI master-worker (MM Borg) |
| 08 | `08_drb_multisite_moea.py` | §7.1–7.2 DRB Case Study | Fig 8 | MPI master-worker (MM Borg) |
| 09 | `09_drb_policy_reeval.py` | §7.3 Policy Re-evaluation | Fig 9 | serial (hand-off to Pywr-DRB) |
| 10 | `10_plot_manuscript_figures.py` | All §/SI figures | Fig 1-9, SI-1, SI-4 | serial; reads every `outputs/expNN_*/` and regenerates PDF figures |

Each `NN_*.slurm` file is a thin wrapper around the matching `NN_*.py`. It sources `_common.sh` for module loads, environment activation, and logging conventions, sets the SLURM topology for the experiment, and calls Python with the appropriate arguments.

## Plotting policy

**All plotting code lives in [`src/plotting/`](../src/plotting/).** Scripts must never define their own plotting functions inline beyond a three-line `fig.savefig(...)` wrapper. Each function in `src/plotting/` is tagged with the manuscript section and figure it produces:

| Module | Figures |
|---|---|
| `src/plotting/style.py` | shared `rcParams`, color palette, water-year month labels |
| `src/plotting/analytic.py` | Fig 1 (concept), Fig 2 (2D/3D analytic), Fig 3 (ε-NFE heatmap), SI-1 (hyperplane) |
| `src/plotting/drought_space.py` | Fig 5 (Kirsch Pareto), Fig 8 (DRB multi-site) |
| `src/plotting/coverage.py` | Fig 7 (library vs MOEA-FIND coverage) |
| `src/plotting/trace_diagnostics.py` | Fig 6 (plausibility: acf, FDC, seasonal cycle, Hurst) |
| `src/plotting/convergence.py` | SI-4 (Borg convergence diagnostics) |

Individual experiment scripts (01-09) call these when invoked with `--plot` and write their working PDFs to `figures/`. Script `10_plot_manuscript_figures.py` is the **single source of truth for the publication figure set**: it reads every available `outputs/expNN_*/` directory and regenerates `figures/fig*.pdf` from scratch. Any figure that appears in the manuscript must be produced by script 10.

## Conventions

- **Outputs.** Each script writes to `outputs/expNN_<slug>/` with `results.json`, any `*.npz` artifacts, and a `config.json` recording the invocation. `--plot` additionally writes a PDF figure to `figures/figNN_<slug>.pdf`.
- **Determinism.** Every driver takes `--seed`. Multi-seed experiments pass `--seeds` as a space-separated list.
- **SynHydro dependency.** Scripts that need SynHydro import it at runtime. Analytic drivers (01, 02, 03) stub it so they can run in a bare Python environment.
- **Borg.** Scripts tagged *MM Borg* expect `borg.py` and compiled binaries on `PYTHONPATH`. See `REPRODUCE.md` (TBD) for the Borg license and build workflow on HPC.

## Local vs HPC

- **Local smoke tests** (workstation): `python scripts/NN_*.py --help` for the CLI, run with small `--nfe` and `--seeds`.
- **HPC production runs**: `sbatch scripts/NN_*.slurm` after editing `_common.sh` for your cluster's module and partition names.

## Cluster setup

`_common.sh` centralizes:
- `#SBATCH` account/partition defaults (override in each `.slurm` as needed).
- Module loads (`module load python mpi`).
- Virtualenv activation.
- `PYTHONPATH` export so `src/` imports resolve.
- `slurm_logs/` directory creation.

Edit the `CLUSTER_*` variables at the top of `_common.sh` once per HPC; individual `.slurm` files should not need per-cluster edits.

## Numbering policy

Scripts are numbered in *manuscript narrative order*, not chronological. If a new experiment is inserted, it gets a `NNa_` suffix rather than renumbering (e.g., `04a_kirsch_single_site_alt_threshold.py`). Section reordering in a WRR methods paper is rare once the outline is fixed; the numeric prefix buys you the "where does this figure come from" lookup at `ls` time.
