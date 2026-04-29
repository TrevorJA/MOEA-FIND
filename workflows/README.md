# MOEA-FIND workflows

Stage-organized experimental pipeline for the MOEA-FIND manuscript.
Each stage is a self-contained folder with its compute drivers, paired
plotting drivers, SLURM batch files, and a `README.md` explaining how
to run that stage. Open the per-stage README to find everything you
need; do not jump between sibling stage folders.

## Stages

| Stage | Manuscript | Purpose | Run order |
|---|---|---|---|
| [01_analytic_validation](01_analytic_validation/) | Figs 1-4, SI-1     | L1 + epsilon-tile proof on synthetic objective space; K=2..6 hypercube dimension sweep | Standalone |
| [02_calibration](02_calibration/)                 | SI-A..E            | Bootstrap tolerances for hydrologic constraints, DV-uniformity, wrapper fidelity / geometry, Kirsch convergence, metric blocks | Before 03, 04 |
| [03_kirsch_library](03_kirsch_library/)           | Fig 7              | 10K-trace Kirsch library + LHS / Sobol / random subsample baselines                | After 02; before 04 |
| [04_moea_find_single_site](04_moea_find_single_site/) | Figs 5-7, SI-F..H | Production single-site Cannonsville MOEA-FIND + ablations (wrapper-mode, DV-uniformity, event-level)            | After 02, 03 |
| [06_pywrdrb_reeval](06_pywrdrb_reeval/)               | Section 7.3, SI-I | Pywr-DRB simulation across the Pareto archive + drought-coverage verification. `replay_pareto_to_multisite_monthly` reuses Kirsch monthly indexes across DRB sites and preserves spatial correlation by construction; the multi-site case study does not need a separate Borg run. | After 04 |
| [07_scenario_discovery](07_scenario_discovery/)       | Fig 9, SI-J        | Satisficing labels + GBT classifiers over hazard-feature space          | After 06 |
| [08_nyc_sensitivity](08_nyc_sensitivity/)             | SI-12              | Delta + PAWN + RBD-FAST sensitivity of NYC outcomes to hazard characteristics | After 06 |
| [99_manuscript_figures](99_manuscript_figures/)       | All                | Single source of truth for `figures/main/` and `figures/supplementary/` | After all upstream stages |

The numeric prefix preserves manuscript-narrative ordering. Stage 05
(multi-site Borg) was removed because the per-site DV vector explosion
adds no information beyond what `policy_reeval.py` already produces.

## Conventions

### Compute / plot split

Compute drivers write only numerical artifacts (`config.json`,
`results.json`, npz, parquet) under
`outputs/<stage>/<driver>/<slug>/`. Figures are produced by paired
plotting drivers under `workflows/<stage>/plots/<driver>.py` that read
those artifacts and write PDFs to
`figures/<stage>/<driver>/<slug>/`. Re-rendering a figure must never
require re-running an experiment.

### SLURM is self-contained

Every `.slurm` file bakes its own primary CLI arguments. Variants are
expressed by maintaining multiple slurm files (e.g.
`run_moea_find.slurm` vs a separate smoke variant) or by SLURM arrays
with a static lookup table inside the slurm. The `--export=ALL,...`
pattern is not used; environment-variable overrides are not part of
the production contract.

Logs land in `workflows/slurm_logs/` (gitignored). `_common.sh` at the
top of `workflows/` is sourced by every slurm via
`${SLURM_SUBMIT_DIR}/workflows/_common.sh`.

### Paths and slugs

- `src/paths.py` resolves output and figure directories:
  - `stage_output_dir(stage, driver, slug=None)` ->
    `outputs/<stage>/<driver>[/<slug>]/`
  - `stage_figure_dir(stage, driver, slug=None)` ->
    `figures/<stage>/<driver>[/<slug>]/`
  - `manuscript_figure_dir(kind="main"|"supplementary")` ->
    `figures/<kind>/` (only stage 99 writes here).
- Slugs are stage-specific helpers under `src/slugs.py` (e.g.
  `moea_slug`, `library_slug`, `subsample_slug`, `analytic_slug`,
  `pywrdrb_slug`). They are deterministic in seeds.
- Determinism: every compute driver takes a `--seed`; the slug encodes
  the seed.

### Parallelism budget

Total concurrent core demand must fit within 10 nodes x 40 cores =
400 cores. Ablation arrays use `%N` throttles to enforce this — see
each stage README for the per-stage table.

### Borg MOEA topology

Every optimization run dispatches MM Borg via MPI (DD-07 -- the
EpsNSGAII platypus stand-in is removed). `src.borg_runner.run_borg_mm`
auto-picks `n_islands` from the MPI rank count using
``max(2, n_ranks // 16)`` for any allocation of >= 5 ranks; single-island
master-slave mode SIGSEGVs on a worker rank during MPI startup once
the per-master worker count exceeds ~6, and using >= 2 islands is the
robust fix. The wrapper also divides the caller's ``nfe`` by
``n_islands`` because Borg's ``solveMPI(maxEvaluations=...)`` is the
per-island budget, not the total. See
``docs/borg_integration_notes.md`` for the full bug analysis.

Concurrent-mpirun stress: when many SLURM array tasks each launch their
own ``mpirun -np N`` simultaneously on the same compute node, Borg MM
elevates its segfault rate. Stage 01 array slurms throttle to 8
concurrent cells (``--array=0-N%8``) to keep failures rare.

## SI subsection mapping

| SI subsection | Driver(s) |
|---|---|
| SI-A Constraint calibration         | `02_calibration/constraint_calibration.py` |
| SI-B DV-uniformity calibration      | `02_calibration/dv_uniformity_calibration.py` |
| SI-C Wrapper fidelity + geometry    | `02_calibration/wrapper_fidelity.py` + `wrapper_geometry.py` |
| SI-D Kirsch convergence             | `02_calibration/kirsch_convergence.py` |
| SI-E Metric-set historical blocks   | `02_calibration/metric_blocks.py` |
| SI-F Wrapper-mode ablation          | `04_moea_find_single_site/wrapper_mode_ablation.py` |
| SI-G DV-uniformity ablation         | `04_moea_find_single_site/dv_uniformity_ablation.py` |
| SI-H Event-level Kirsch             | `04_moea_find_single_site/event_level.py` |
| SI-I Drought-coverage verification  | `06_pywrdrb_reeval/verify_drought_coverage.py` |
| SI-J Satisficing manifold + GBT     | `07_scenario_discovery/satisficing_sweep.py` + `scenario_discovery_plots.py` |
| SI-12 NYC sensitivity               | `08_nyc_sensitivity/run_sa.py` + `compare_methods.py` |

## Local vs HPC

- **Local smoke test:** `python workflows/<stage>/<driver>.py --help`
  to verify the CLI parses; run analytic stages with small `--nfe`
  (these work in a bare Python environment without SynHydro).
- **HPC production:** `sbatch workflows/<stage>/slurm/<driver>.slurm`.
  Edit `_common.sh` for your cluster's module and partition names.
