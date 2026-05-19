# MOEA-FIND workflows

Stage-organized experimental pipeline for the MOEA-FIND manuscript.
Each numbered stage is a self-contained **compute** step: its driver
scripts write only numerical artifacts. All **figure** scripts live in
the two stage-99 figure folders, not inside the per-stage folders.

## Stages

| Stage | Manuscript | Purpose | Run order |
|---|---|---|---|
| [01_analytic_validation](01_analytic_validation/) | Fig 4, SI-1 | L1 + epsilon-tile proof on synthetic objective space; K=2..6 dimension sweep; eps×NFE sweep | Standalone |
| [02_calibration](02_calibration/) | SI-A/B + DD-15 (K*,T*) | Bootstrap tolerances for hydrologic & DV-uniformity constraints; T-sensitivity / decision-matrix; short-block & bounded-metric selection; first-event envelope | Before 03, 04 |
| [03_kirsch_library](03_kirsch_library/) | Fig 7 | 10K-trace Kirsch library + LHS/Sobol/random subsample baselines | After 02; before 04 |
| [04_moea_find_single_site](04_moea_find_single_site/) | Figs 5-7, SI-F/G/H | Production single-site Cannonsville MOEA-FIND + ablations | After 02, 03 |
| [06_pywrdrb_reeval](06_pywrdrb_reeval/) | §3.3, SI-I | Pywr-DRB simulation across the Pareto archive + coverage verification | After 04 |
| [07_scenario_discovery](07_scenario_discovery/) | Fig 7, SI-J | Satisficing labels + GBT classifiers over hazard-feature space | After 06 |
| [08_nyc_sensitivity](08_nyc_sensitivity/) | SI-12 | Delta + PAWN + RBD-FAST sensitivity of NYC outcomes | After 06 |
| [09_magnitude_varying_sa](09_magnitude_varying_sa/) | SI-12b | Per-percentile magnitude-varying sensitivity | After 06 |
| [99_manuscript_figures](99_manuscript_figures/) | Main text | The single job that builds every **main-text** figure (`figures/main/`) | After all upstream |
| [99_supporting_info](99_supporting_info/) | SI-C/D | SI-only **compute** diagnostics (Kirsch convergence, wrapper geometry/fidelity) | After 03/04 |
| [99_supporting_info_figures](99_supporting_info_figures/) | SI | The single job that builds every **SI** figure (`figures/supplementary/`) | After upstream |

Stage 05 (multi-site Borg) was removed: the per-site DV explosion adds
nothing beyond `policy_reeval.py`. The numeric prefix preserves
manuscript-narrative order.

## Conventions

### Compute vs plotting (post spring-clean reorg)

- **Reusable plotting *functions*** live in `src/plotting/*.py` (library
  modules — `drought_space`, `sensitivity`, `trace_diagnostics`, …).
- **Scripts that *call* those functions** to render a figure live only
  in `workflows/99_manuscript_figures/` (main-text figures) or
  `workflows/99_supporting_info_figures/` (SI figures). There are **no**
  per-stage `plots/` folders and **no** per-stage plotting SLURM jobs.
- Compute drivers write `outputs/<stage>/<driver>/<slug>/` only;
  re-rendering a figure never re-runs an experiment.

### SLURM: ≤2 per step, no variant sprawl

- Each step has **one** primary `.slurm` (plus at most a `*_smoke`
  second). Old `_v2/_coarse/_fixed/_iter*/_ks/_index/_T5` variant slurms
  were removed.
- Framings are expressed as **config files**, not as duplicated slurm.
  Stage 04 is fully parameterised:
  `sbatch workflows/04_moea_find_single_site/slurm/run_moea_find.slurm [CONFIG] [SEED] [SUFFIX]`
  where `CONFIG` is a YAML basename under
  `workflows/04_moea_find_single_site/configs/`. Every prior framing
  (short_block_drb_v2, first_event_ssi3_t10, primary, …) remains
  runnable through that one slurm.
- All figures are built by exactly two jobs:
  `99_manuscript_figures/slurm/make_figures.slurm` and
  `99_supporting_info_figures/slurm/make_si_figures.slurm`.
- `_common.sh` (top of `workflows/`) is sourced by every slurm; logs go
  to `workflows/slurm_logs/` (gitignored).

### Paths and slugs

`src/io_paths/paths.py` resolves directories:
- `stage_output_dir(stage, driver, slug=None)` → `outputs/<stage>/<driver>[/<slug>]/`
- `stage_figure_dir(stage, driver, slug=None)` → `figures/<stage>/<driver>[/<slug>]/`
- `manuscript_figure_dir(kind="main"|"supplementary")` → `figures/<kind>/`

Slugs are the deterministic, seed-encoding helpers in
`src/io_paths/slugs.py` (`moea_slug`, `library_slug`, `subsample_slug`,
`analytic_slug`, `pywrdrb_slug`). Stage 04 builds its output slug via
`moea_slug` (the legacy `make_variant_slug` was retired).

### Borg MOEA topology

Every optimization run dispatches MM Borg via MPI. `run_moea_find.slurm`
is `--exclusive` at full rank; smokes drop NFE, not ranks. See
`docs/borg_integration_notes.md`. Concurrent-mpirun stress: stage-01
array slurms throttle (`--array=0-N%8`).

## SI subsection mapping

| SI subsection | Compute driver | Figure script (in 99_*_figures) |
|---|---|---|
| SI-A Constraint calibration   | `02_calibration/constraint_calibration.py` | `constraint_calibration.py` |
| SI-B DV-uniformity calibration| `02_calibration/dv_uniformity_calibration.py` | `dv_uniformity_calibration.py` |
| SI-C Wrapper fidelity+geometry| `99_supporting_info/wrapper_fidelity.py` + `wrapper_geometry.py` | `wrapper_fidelity.py` + `wrapper_geometry.py` |
| SI-D Kirsch convergence       | `99_supporting_info/kirsch_convergence.py` | `kirsch_convergence.py` |
| SI-F Wrapper-mode ablation    | `04_moea_find_single_site/wrapper_mode_ablation.py` | `wrapper_mode_ablation.py` |
| SI-G DV-uniformity ablation   | `04_moea_find_single_site/dv_uniformity_ablation.py` | `dv_uniformity_ablation.py` |
| SI-H Event-level Kirsch       | `04_moea_find_single_site/event_level.py` | `first_event_archive_scatter.py` |
| SI-I Drought-coverage verify  | `06_pywrdrb_reeval/verify_drought_coverage.py` | `verify_drought_coverage.py` |
| SI-J Satisficing + GBT        | `07_scenario_discovery/satisficing_sweep.py` | `scenario_discovery_plots.py`* |
| SI-12 NYC sensitivity         | `08_nyc_sensitivity/run_sa.py` | `run_sa.py`, `compare_methods.py` |
| SI-12b Magnitude-varying SA   | `09_magnitude_varying_sa/run_mv_sa.py` | `run_mv_sa.py` |

\* `scenario_discovery_plots.py` also feeds main-text Fig 7, so it lives
in `99_manuscript_figures/`. The former SI-E metric-block exploration
(`metric_explorer`, `metric_blocks`) was deleted as redundant —
superseded by the DD-15 `t_sensitivity_*` pipeline. The SI is not yet
final; more subsections may be added.

## Local vs HPC

- **Local smoke:** `python workflows/<stage>/<driver>.py --help`
  (analytic stages run in a bare env without SynHydro).
- **HPC:** `sbatch workflows/<stage>/slurm/<driver>.slurm`. Edit
  `_common.sh` for your cluster's modules/partition.
