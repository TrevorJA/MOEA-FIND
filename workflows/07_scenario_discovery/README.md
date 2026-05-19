# Stage 07 -- Scenario discovery


## Purpose

Post-process the per-realization metric bank from stage 06 into binary
satisficing labels and gradient-boosted classifiers over drought-feature
space. Cheap to re-run -- adding or tuning a satisficing definition does
not require re-running Pywr-DRB. Backs the manuscript's main-text
scenario-discovery figure (Fig 9) and the SI-J classifier diagnostics.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [satisficing_sweep.py](satisficing_sweep.py)                                 | Per-definition GBT classifier sweep over [satisficing_manifest.yaml](satisficing_manifest.yaml). | SI-J  |
| [plots/scenario_discovery_plots.py](plots/scenario_discovery_plots.py)       | Final SD figure generator (single + multi-panel maps, AUC bars, gbt + logreg). | Fig 9 |

`scenario_discovery_plots.py` is treated as the plotting driver: it
re-fits classifiers and refreshes the satisficing table on demand
(cheap), but writes those numerical derivatives under its own outputs
slug rather than mutating stage 06.

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/07_scenario_discovery/<driver>/<slug>/`. Figures land under
`figures/07_scenario_discovery/<driver>/<slug>/`. Re-rendering a
figure never requires re-running stage 06.

## Run order

Must run **after** stages 04 (Pareto archive) and 06 (Pywr-DRB
re-evaluation -> metric bank). Both drivers consume:

- `outputs/04_moea_find_single_site/run_moea_find/<slug>/results.json`
  (Pareto archive + drought characteristics)
- `outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet`
  (and `drought_levels.json` alongside it)

By default `<src_slug>` is the directory name of the Pareto archive's
parent; override with `--src-slug` when stage 06 was re-evaluated under
a different slug.

## Slurm

| Compute slurm | Plotting slurm |
|---|---|
| `slurm/satisficing_sweep.slurm`   | `slurm/scenario_discovery_plots.slurm` |

All slurm scripts are self-contained -- arguments are baked into the
slurm file. Variant changes are made by editing the slurm script, not
via `--export`.

## Outputs

```
outputs/07_scenario_discovery/
  satisficing_sweep/<slug>/
    {moea_find,baseline?}/{classifiers/<def>/, figures/boundary_<def>.pdf}
    classifier_summary.csv
    figures/manifest_summary.pdf
  scenario_discovery_plots/<slug>/
    classification.csv
    drought_levels.json
    {gbt,logreg}/{classifiers/<def>/, classifier_summary.csv}
```

## Figures

```
figures/07_scenario_discovery/
  scenario_discovery_plots/<slug>/
    fig09_satisficing_map.pdf
    fig09_satisficing_map_{gbt,logreg}.pdf
    fig_manifest_summary_{gbt,logreg}.pdf
```

Final manuscript PDFs are regenerated from these working figures by
[workflows/99_manuscript_figures/make_figures.py](../99_manuscript_figures/make_figures.py).
