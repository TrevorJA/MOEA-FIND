# 99_supporting_info

SI-only **compute** diagnostics relocated out of `02_calibration` during
the spring-cleaning reorg, because they are not on the production
critical path (nothing downstream consumes their outputs) but remain
valuable small tests for the manuscript SI.

| Driver | SI | Purpose |
|---|---|---|
| `kirsch_convergence.py` | SI-D | Kirsch ensemble convergence vs ensemble size × seed |
| `wrapper_geometry.py`   | SI-C | Single-DV sweep / seasonal-profile geometry diagnostic |
| `wrapper_fidelity.py`   | SI-C | Ensemble + drought-cloud fidelity vs history/Kirsch |

## Run

One job runs all three (compute only; figures come from the matching
scripts in `workflows/99_supporting_info_figures/`):

```
sbatch workflows/99_supporting_info/slurm/supporting_info.slurm [LIBRARY_NPZ] [MOEA_FRONT_JSON]
```

The two optional args point `kirsch_convergence` at a stage-03 library
`characteristics.npz` and a stage-04 `results.json`; update them to
match the slugs produced by your reruns.
