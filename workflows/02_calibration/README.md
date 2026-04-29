# Stage 02 -- Calibration

## Purpose

Bootstrap-calibrated tolerance files and diagnostic figures consumed by
stages 03 and 04. Each driver here backs one SI subsection.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [constraint_calibration.py](constraint_calibration.py)         | Bootstrap tolerances for the hydrologic-constraint mode.       | SI-A |
| [dv_uniformity_calibration.py](dv_uniformity_calibration.py)   | Bootstrap tolerance for the DV-uniformity constraint (residual + index). | SI-B |
| [wrapper_fidelity.py](wrapper_fidelity.py)                     | Wrapper preserves Kirsch marginals + drought-space coverage under U[0,1] DVs. | SI-C |
| [wrapper_geometry.py](wrapper_geometry.py)                     | Mapping smoothness sweep; index vs residual cartoon.            | SI-C |
| [kirsch_convergence.py](kirsch_convergence.py)                 | Convergence of drought-characteristic range / coverage with library size. | SI-D |
| [metric_blocks.py](metric_blocks.py)                           | Per-metric variability across historical T-year blocks.        | SI-E |

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/02_calibration/<driver>/[<slug>/]`. Figures are produced by
the paired plotting drivers under [plots/](plots/) which read those
artifacts and write to `figures/02_calibration/<driver>/`. Re-rendering
a figure never requires re-running the bootstrap.

## Slurm

| Compute slurm | Plotting slurm |
|---|---|
| `slurm/constraint_calibration.slurm`     | `slurm/plots/constraint_calibration.slurm`     |
| `slurm/dv_uniformity_calibration.slurm`  | `slurm/plots/dv_uniformity_calibration.slurm`  |
| `slurm/wrapper_fidelity.slurm`           | `slurm/plots/wrapper_fidelity.slurm`           |
| `slurm/wrapper_geometry.slurm`           | `slurm/plots/wrapper_geometry.slurm`           |
| `slurm/kirsch_convergence.slurm`         | `slurm/plots/kirsch_convergence.slurm`         |
| `slurm/metric_blocks.slurm`              | `slurm/plots/metric_blocks.slurm`              |

`dv_uniformity_calibration` is a 2-task array (residual + index modes).
All other compute slurms are single jobs (1-16 cores).

All slurms are self-contained -- arguments are baked into the slurm
files, no `--export=ALL` overrides.

## Outputs

```
outputs/02_calibration/
  constraint_calibration/{config.json, calibrated_tolerances.json, ...}
  dv_uniformity_calibration/{residual,index}/{config.json, calibrated_dv_tolerances.json, ...}
  wrapper_fidelity/{config.json, ensembles_2d.npz, drought_clouds.npz, ...}
  wrapper_geometry/{config.json, sweep.npz, geometry_summary.json}
  kirsch_convergence/{config.json, convergence.json}
  metric_blocks/{config.json, block_chars.csv, per_metric_summary.csv,
                 per_preset_summary.csv, full_hist_chars.pkl}
```

## Figures

```
figures/02_calibration/
  constraint_calibration/...
  dv_uniformity_calibration/{residual,index}/...
  wrapper_fidelity/...
  wrapper_geometry/...
  kirsch_convergence/...
  metric_blocks/...
```

## Run order

Stage 02 must complete before stage 04 ablation runs. The two JSON
tolerance files consumed downstream are:
- `outputs/02_calibration/constraint_calibration/calibrated_tolerances.json`
- `outputs/02_calibration/dv_uniformity_calibration/{residual,index}/calibrated_dv_tolerances.json`
