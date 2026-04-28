# MOEA-FIND workflows

Stage-organized experimental pipeline for the MOEA-FIND manuscript.
Each stage is a self-contained folder with its driver scripts, YAML
configs, SLURM batch files, and a `README.md` that explains how to run
that stage. Open the per-stage README to find everything you need —
no jumping between sibling folders.

## Stages

| Stage | Section(s) | Purpose | Run order |
|---|---|---|---|
| [01_analytic_validation](01_analytic_validation/) | §5, SI-1 | Manhattan-norm proof on synthetic objective space; K=1..6 dimension sweep | Standalone |
| [02_calibration](02_calibration/)                 | SI-A..E   | Constraint, DV-uniformity, wrapper, Kirsch convergence calibrations    | Before 03, 04 |
| [03_kirsch_library](03_kirsch_library/)           | §6.3 (Fig 7) | Library + LHS/Sobol subsample baseline                              | After 02, before 04 |
| [04_moea_find_single_site](04_moea_find_single_site/) | §6.1-6.4, SI-F..H | Core single-site MOEA-FIND + ablations                            | After 02, 03 |
| [05_moea_find_multisite](05_moea_find_multisite/)     | §7.1-7.2 (Fig 8) | DRB multi-site application                                       | After 02 |
| [06_pywrdrb_reeval](06_pywrdrb_reeval/)               | §7.3, SI-I    | Pywr-DRB simulation on Pareto droughts (simulation only)         | After 04 or 05 |
| [07_scenario_discovery](07_scenario_discovery/)       | §7.3 (Fig 9), SI-J | Satisficing labels + GBT classifiers (post-processing)      | After 06 |
| [08_nyc_sensitivity](08_nyc_sensitivity/)             | (planned)     | Exploratory-modeling sensitivity of NYC reservoir outcomes       | Placeholder; later session |
| [99_manuscript_figures](99_manuscript_figures/)       | All           | Regenerate `figures/main/` and `figures/supplementary/` from upstream outputs | After all upstream stages |

The numeric prefix preserves manuscript-narrative ordering. New stages
inserted into the sequence get a `0Na_` suffix rather than renumbering
(per the policy from the previous workflow layout).

## Conventions

- **Drivers.** One Python `.py` file per experiment, located inside
  the stage folder (not in a separate `experiments/` directory).
- **SLURM.** Each stage has a `slurm/` subfolder with one `.slurm`
  per driver. Cluster-debugging SLURM files without a matching `.py`
  live under `workflows/_scratch/slurm/` and are not part of the
  reproducible pipeline. `_common.sh` (module loads, venv activation,
  log directory) stays at the top level and is sourced by every
  SLURM script via `${SLURM_SUBMIT_DIR}/workflows/_common.sh`.
- **Logs.** All SLURM stdout/stderr land in `workflows/slurm_logs/`
  (gitignored). One central directory across all stages.
- **YAML configs.** Stage 04 supports YAML presets for ablation arms
  (`workflows/04_moea_find_single_site/configs/*.yaml`). The driver
  takes `--config <path>`; CLI flags still override loaded values.
- **Outputs.** Every driver writes to `outputs/expNN_<slug>/<variant>/`
  with `config.json` (invocation record), `results.json` (summary),
  and stage-specific artifacts. See [outputs/README.md](../outputs/README.md)
  for the slug catalog.
- **Slugs.** Variant slugs are built by `src.slugs.*` (preferred) or
  the legacy `src.experiment_utils.make_variant_slug` (kept for
  back-compat with existing output directories). Format:
  `stage__key=val__key=val__s=seed`. See `src/slugs.py` for the
  stage-specific helpers.
- **Determinism.** Every driver takes `--seed`; output slugs encode it.
- **Figures.** Working figures land in `figures/<stage>/`. Confirmed
  manuscript figures live in `figures/main/` and
  `figures/supplementary/` and are regenerated only by stage 99.

## Plotting policy

All plotting code lives in [`src/plotting/`](../src/plotting/). Drivers
must not define inline plotting beyond a thin `fig.savefig(...)` wrapper.
See `99_manuscript_figures/README.md` for the figure promotion contract.

## SI subsection mapping

Each method-consideration script earns SI text. The mapping below is
mirrored in [manuscript/drafts/](../manuscript/drafts/) as one stub per
subsection.

| SI subsection | Driver |
|---|---|
| §SI-A Constraint calibration         | `02_calibration/constraint_calibration.py` |
| §SI-B DV-uniformity calibration      | `02_calibration/dv_uniformity_calibration.py` |
| §SI-C Wrapper fidelity & geometry    | `02_calibration/wrapper_fidelity.py` + `wrapper_geometry.py` |
| §SI-D Kirsch convergence             | `02_calibration/kirsch_convergence.py` |
| §SI-E Metric-set historical blocks   | `02_calibration/metric_blocks.py` |
| §SI-F Wrapper-mode ablation          | `04_moea_find_single_site/wrapper_mode_ablation.py` |
| §SI-G DV-uniformity ablation         | `04_moea_find_single_site/dv_uniformity_ablation.py` |
| §SI-H Event-level Kirsch             | `04_moea_find_single_site/event_level.py` |
| §SI-I Drought-coverage verification  | `06_pywrdrb_reeval/verify_drought_coverage.py` |
| §SI-J Satisficing manifold + GBT     | `07_scenario_discovery/satisficing_sweep.py` + `scenario_discovery_plots.py` |

Each subsection targets ~2-3 multi-panel SI figures (clean, multi-panel
academic style). The secondary scoping task tightens these to budget.

## Local vs HPC

- **Local smoke test:** `python workflows/0N_<stage>/<driver>.py --help`
  for the CLI; run with small `--nfe` (analytic stages run on a bare
  Python environment without SynHydro).
- **HPC production:**
  `sbatch workflows/0N_<stage>/slurm/<driver>.slurm` after editing
  `_common.sh` for your cluster's module and partition names.
