# Stage 99 - Manuscript figure assembly

## Purpose

Single source of truth for the final manuscript figure set. This is the
**only** workflow allowed to write into `figures/main/` and
`figures/supplementary/`. All other workflows write per-stage figures
under `figures/<stage>/<driver>/<slug>/` via `src.paths.stage_figure_dir`;
stage 99 promotes a curated subset into the manuscript folders.

## Inputs

Reads from the restructured outputs tree at `outputs/<stage>/<driver>/<slug>/`:

| Figure | Source |
|---|---|
| fig04 | `outputs/01_analytic_validation/dimension_sweep/` |
| fig05, fig06 | `outputs/04_moea_find_single_site/run_moea_find/<slug>/` |
| fig06 (Kirsch library) | `outputs/03_kirsch_library/build_library/` |
| fig07 | `outputs/06_pywrdrb_reeval/policy_reeval/<slug>/` |
| figSI01 | `outputs/01_analytic_validation/analytic_2d/`, `.../analytic_3d/` |

fig01, fig02, fig03 carry no upstream data dependency (placeholder /
schematic / Manhattan construction figure).

## Outputs

- `figures/main/figNN_<descriptor>.pdf`
- `figures/supplementary/figSI<NN>_<descriptor>.pdf`

Resolved through `src.paths.manuscript_figure_dir("main")` and
`src.paths.manuscript_figure_dir("supplementary")`.

## CLI

```bash
# Assemble every figure (fig05/06/07 SKIP without --single-site-slug).
python workflows/99_manuscript_figures/make_figures.py

# Promote a specific MOEA-FIND single-site variant.
python workflows/99_manuscript_figures/make_figures.py \
    --single-site-slug residual_T20_nfe200000_s42_constrained

# Regenerate one figure.
python workflows/99_manuscript_figures/make_figures.py --only fig04
```

If a contributing upstream output is missing, the corresponding figure
is skipped with a warning rather than failing the run.

## SLURM

`sbatch workflows/99_manuscript_figures/slurm/make_figures.slurm`

The slug fed to `--single-site-slug` is baked into the SLURM script;
edit it there when promoting a different variant to the manuscript.
