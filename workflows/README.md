# MOEA-FIND workflows

Stage-organized experimental pipeline for the MOEA-FIND manuscript.
Each numbered stage is a self-contained **compute** step that writes
numerical artifacts only; figure rendering lives in the two stage-99
folders.

## Stages

| Stage | Manuscript | Purpose | Run order |
|---|---|---|---|
| [01_analytic_validation](01_analytic_validation/) | Fig 4, SI-1 | K-dim analytic benchmark sweep (K=2..6); ε × NFE sensitivity grid | Standalone |
| [02_calibration](02_calibration/) | production setup | Anderson-Darling DV-uniformity threshold + first-event envelope diagnostic | Before 03, 04 |
| [03_kirsch_library](03_kirsch_library/) | Fig 7 | Kirsch library generation → metric characterization → LHS/Sobol/random subsample baselines | After 02; before 04 |
| [04_moea_find_single_site](04_moea_find_single_site/) | Figs 5-7, SI-F/H | Production single-site Cannonsville MOEA-FIND (first-event T=10y default; primary + wrapper ablation configs) | After 02, 03 |
| [06_pywrdrb_reeval](06_pywrdrb_reeval/) | §3.3, SI-I | Pywr-DRB simulation across the Pareto archive (prepare → simulate → aggregate) + coverage verification | After 04 |
| [07_scenario_discovery](07_scenario_discovery/) | Fig 7, SI-J | Satisficing labels + GBT classifiers over the hazard-feature space | After 06 |
| [08_nyc_sensitivity](08_nyc_sensitivity/) | SI-12 | Delta / PAWN / RBD-FAST sensitivity of NYC outcomes | After 06 |
| [09_magnitude_varying_sa](09_magnitude_varying_sa/) | SI-12b | Per-percentile magnitude-varying sensitivity | After 06 |
| [99_manuscript_figures](99_manuscript_figures/) | Main text | Single job that builds every main-text figure (`figures/main/`) | After all upstream |
| [99_supporting_info](99_supporting_info/) | SI-C/D | SI-only compute diagnostics (Kirsch convergence, wrapper geometry / fidelity) | After 03/04 |
| [99_supporting_info_figures](99_supporting_info_figures/) | SI | Single job that builds every SI figure (`figures/supplementary/`) | After upstream |

The numeric prefix preserves manuscript-narrative order; stage 05
(multi-site Borg) is intentionally absent.

## Conventions

### Compute vs plotting

- Reusable plotting **functions** live in `src/plotting/*.py`.
- Scripts that **call** those functions to render a figure live only
  in `99_manuscript_figures/` (main text) or
  `99_supporting_info_figures/` (SI). There are no per-stage plotting
  scripts and no per-stage figure SLURMs.
- Compute drivers write `outputs/<stage>/<driver>/<slug>/`; re-rendering
  a figure never re-runs an experiment.

### SLURM: one primary + at most one smoke

- Each step has **one** primary `.slurm` (plus optionally a `*_smoke`).
- Framings are expressed as YAML configs, not duplicated SLURM files.
  Stage 04 is fully parameterised:
  `sbatch workflows/04_moea_find_single_site/slurm/run_moea_find.slurm [CONFIG] [SEED] [SUFFIX]`
  where `CONFIG` is a YAML basename under
  `workflows/04_moea_find_single_site/configs/`.
- All figures are built by exactly two jobs:
  `99_manuscript_figures/slurm/make_figures.slurm` and
  `99_supporting_info_figures/slurm/make_si_figures.slurm`.
- `_common.sh` (top of `workflows/`) is sourced by every SLURM; logs go
  to `workflows/slurm_logs/` (gitignored).

### Paths and slugs

`src/io_paths/paths.py` resolves directories:
- `stage_output_dir(stage, driver, slug=None)` → `outputs/<stage>/<driver>[/<slug>]/`
- `stage_figure_dir(stage, driver, slug=None)` → `figures/<stage>/<driver>[/<slug>]/`
- `manuscript_figure_dir(kind="main"|"supplementary")` → `figures/<kind>/`

Slugs come from the canonical helpers in `src/io_paths/slugs.py`:
`moea_slug`, `library_slug`, `subsample_slug`, `analytic_slug`,
`pywrdrb_slug`. Hand-rolled slug strings are not used by drivers.

### Borg MOEA topology

Every optimization run dispatches MM Borg via MPI.
`run_moea_find.slurm` is `--exclusive` at full rank; smokes drop NFE,
not ranks. See `docs/borg_integration_notes.md`. Concurrent-mpirun
stress: stage-01 array slurms throttle (`--array=0-N%8`).

## SI subsection mapping

| SI subsection | Compute driver | Figure script (in 99_supporting_info_figures/) |
|---|---|---|
| SI-C Wrapper fidelity + geometry | `99_supporting_info/wrapper_fidelity.py` + `wrapper_geometry.py` | `wrapper_fidelity.py` + `wrapper_geometry.py` |
| SI-D Kirsch convergence          | `99_supporting_info/kirsch_convergence.py` | `kirsch_convergence.py` |
| SI-F Wrapper-mode ablation       | `04_moea_find_single_site/wrapper_mode_ablation.py` | `wrapper_mode_ablation.py` |
| SI-H First-event archive         | `04_moea_find_single_site/run_moea_find.py --config first_event_ssi3_t10.yaml` | `first_event_archive_scatter.py` |
| SI-I Drought-coverage verify     | `06_pywrdrb_reeval/verify_drought_coverage.py` | `verify_drought_coverage.py` |
| SI-J Satisficing + GBT           | `07_scenario_discovery/satisficing_sweep.py` | `scenario_discovery_plots.py`* |
| SI-12 NYC sensitivity            | `08_nyc_sensitivity/run_sa.py` | `run_sa.py`, `compare_methods.py` |
| SI-12b Magnitude-varying SA      | `09_magnitude_varying_sa/run_mv_sa.py` | `run_mv_sa.py` |

\* `scenario_discovery_plots.py` also feeds main-text Fig 7, so it
lives in `99_manuscript_figures/`.

## Stage 06 prep / sim / aggregate split

`06_pywrdrb_reeval` runs as three independent jobs so the heavy
Kirsch-replay / Nowak-disagg / KDE chain caches on disk and the
Pywr-DRB simulation can be re-run cheaply:

1. `sbatch prepare_windowed.slurm [PARETO_RESULTS]` — DV → daily
   multisite HDF5 + per-scenario `window.json` (sim_start, sim_end,
   first-event metadata). Writes
   `scenarios/<sid>/{pywrdrb_inputs/, window.json, status_prep.json}`.
   Idempotent (`--overwrite` to force).
2. `sbatch policy_reeval_windowed.slurm [PARETO_RESULTS]` — reads the
   prepped inputs + window, runs Pywr-DRB over the window, writes
   `scenarios/<sid>/simulations/batch_0.hdf5` + `status_sim.json`.
   Also idempotent.
3. `sbatch aggregate_windowed.slurm` — concatenates per-scenario
   outputs into `results/metric_bank.parquet` for stages 07/08/09.

## Local vs HPC

- **Local smoke:** `python workflows/<stage>/<driver>.py --help`
  (analytic stages run in a bare env without SynHydro).
- **HPC:** `sbatch workflows/<stage>/slurm/<driver>.slurm`. Edit
  `_common.sh` for your cluster's modules / partition.
