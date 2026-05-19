# Stage 04 -- Single-site MOEA-FIND


## Purpose

Production single-site Cannonsville MOEA-FIND plus the three SI
ablations (wrapper-mode, DV-uniformity, event-level), the post-hoc
library-vs-MOEA baseline comparison, and the post-hoc compare drivers.
Backs the manuscript main results plus SI-F, SI-G, SI-H.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [run_moea_find.py](run_moea_find.py)               | Production single-site MOEA-FIND with MM Borg, residual wrapper, DV-uniformity (AD) constraint. | Figs 5, 6 |
| [baseline_comparison.py](baseline_comparison.py)   | Stage 03 Kirsch library vs stage 04 Pareto coverage. | Fig 7 |
| [dv_uniformity_ablation.py](dv_uniformity_ablation.py) | Hydrologic vs DV-uniformity constraint regimes (per-arm Pareto). | SI-G |
| [wrapper_mode_ablation.py](wrapper_mode_ablation.py)   | Index vs residual wrapper-mode (per-mode Pareto). | SI-F |
| [event_level.py](event_level.py)                       | Short-trace event-level Kirsch (scaffold; --dry-run today). | SI-H |

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/04_moea_find_single_site/<driver>/[<arm-or-mode>/]<slug>/`.
Figures are produced by paired plotting drivers under [plots/](plots/):

| Plot driver | Reads | Writes |
|---|---|---|
| `plots/run_moea_find.py`         | run_moea_find/<slug>/results.json + historical_block_chars.npz | figures/.../run_moea_find/<slug>/ |
| `plots/baseline_comparison.py`   | baseline_comparison/{pooled.npz, comparison_summary.json}      | figures/.../baseline_comparison/  |
| `plots/dv_uniformity_ablation.py`| dv_uniformity_ablation/<arm>/<slug>/results.json               | figures/.../dv_uniformity_ablation/<arm>/<slug>/ |
| `plots/dv_uniformity_compare.py` | every dv_uniformity_ablation/<arm>/<slug>/results.json (SI-G)  | figures/.../dv_uniformity_compare/ |
| `plots/wrapper_mode_ablation.py` | wrapper_mode_ablation/<mode>/<slug>/results.json               | figures/.../wrapper_mode_ablation/<mode>/<slug>/ |
| `plots/wrapper_mode_compare.py`  | every wrapper_mode_ablation/<mode>/<slug>/results.json (SI-F)  | figures/.../wrapper_mode_compare/  |

## Slurm

| Compute slurm | Cores | Plot slurm |
|---|---|---|
| `slurm/event_level.slurm`               | 32  | -- |

All slurms are self-contained: arguments are baked in. Variant changes
are made by editing the slurm or adding a sibling slurm, not by env-var
overrides.

## Run order

1. Stage 02 calibrations must complete first (constraint + DV-uniformity
   for both residual and index modes).
2. Stage 03 `build_library` must complete before `baseline_comparison`.
3. Submit `run_moea_find` for the production main-text result.
4. Optional: submit ablation arrays + their compare plotting jobs.

## Configs

YAML presets under [configs/](configs/) capture canonical argument
bundles (primary, dv_uniformity, wrapper_index, wrapper_residual). The
slurms bake the same argument values directly so the configs are
documentation, not runtime dependencies.

## Outputs

```
outputs/04_moea_find_single_site/
  run_moea_find/<slug>/{config.json, results.json, pareto.npz,
                        historical_block_chars.npz, historical_blocks.npz}
  baseline_comparison/{config.json, comparison_summary.json, pooled.npz}
  dv_uniformity_ablation/<arm>/<slug>/{config.json, results.json, pareto.npz}
  dv_uniformity_compare/{comparison_summary.json}
  wrapper_mode_ablation/<mode>/<slug>/{config.json, results.json, pareto.npz}
  wrapper_mode_compare/{wrapper_comparison_summary.json}
  event_level/<slug>/{config.json}
```

## Figures

```
figures/04_moea_find_single_site/
  run_moea_find/<slug>/{fig05_drought_space_2d.pdf, fig06_drought_space_3d.pdf}
  baseline_comparison/{fig07_library_vs_moea.pdf, fig07_coverage_bars.pdf}
  dv_uniformity_ablation/<arm>/<slug>/drought_space.pdf
  dv_uniformity_compare/{figSI_ablation_*.pdf}
  wrapper_mode_ablation/<mode>/<slug>/drought_space.pdf
  wrapper_mode_compare/{figSI_wrapper_*.pdf}
```
