# Stage 03 — Kirsch library + subsample baseline

## Purpose

Build a large Kirsch library (~10K traces) and produce LHS/Sobol
subsamples in drought-characteristic space as the baseline against which
MOEA-FIND coverage is compared. Backs §6.3 of the manuscript and Fig 7
(library vs MOEA-FIND coverage).

## Run order

Must run **after** stage 02 (uses no calibration output directly, but
uses the same Kirsch wrapper that stage 02 validates). Must run
**before** stage 04's `baseline_comparison.py` driver, which consumes
`outputs/exp06_library_subsample/`.

## Drivers

| Driver | SLURM | Outputs | Working figures |
|---|---|---|---|
| [build_library.py](build_library.py)         | [slurm/build_library.slurm](slurm/build_library.slurm)         | `outputs/exp05_kirsch_library/`    | `figures/kirsch_library/` |
| [subsample_baseline.py](subsample_baseline.py) | [slurm/subsample_baseline.slurm](slurm/subsample_baseline.slurm) | `outputs/exp06_library_subsample/` | `figures/kirsch_library/` |

`build_library.py` consumes a fitted KirschGenerator
(via `src.kirsch_utils.build_kirsch_generator`) and writes
`library.npy` + `characteristics.json` per slug. `subsample_baseline.py`
applies LHS / Sobol selection in the metric-bank space (`src.analysis`).

## Manuscript figures

- Main: contributes data to Fig 7 (library vs MOEA-FIND coverage),
  produced by `99_manuscript_figures/make_figures.py`.
- SI: none planned.

## SI figure budget

No SI section dedicated to this stage; library+subsample is shown only
as the comparison baseline in §6.3.
