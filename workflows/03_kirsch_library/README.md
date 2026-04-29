# Stage 03 -- Kirsch library + subsample baseline

## Purpose

Build a large Kirsch-Nowak library (~10K traces) and produce LHS / Sobol /
random subsamples in drought-characteristic space. The library is the
natural sampling distribution of the Kirsch bootstrap generator; the
subsamples are the baselines against which MOEA-FIND coverage is
compared. Backs the manuscript's main Fig 7
(library vs MOEA-FIND coverage). The actual comparison figure is
rendered by stage 04's `baseline_comparison.py` driver -- stage 03 only
produces the numerical artifacts it consumes.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [build_library.py](build_library.py)         | Generate 10K-trace Kirsch library + per-trace SSI drought metrics. | Fig 7 (data) |
| [subsample_baseline.py](subsample_baseline.py) | LHS / Sobol / random subsample of the library via nearest-neighbor matching. | Fig 7 (data) |

## Compute / plot split

Both drivers are compute-only -- they write numerical artifacts under
`outputs/03_kirsch_library/<driver>/<slug>/` and produce no figures.
The Fig 7 comparison plot is rendered by stage 04's
`baseline_comparison` driver, which reads stage 03 outputs alongside
the MOEA-FIND results.

## Slurm

| Compute slurm |
|---|
| `slurm/build_library.slurm`      |
| `slurm/subsample_baseline.slurm` |

`build_library.slurm` runs serially -- SynHydro handles the ensemble
loop internally. `subsample_baseline.slurm` loops over `{lhs, sobol,
random}` inside one job so all three method outputs are produced
together.

All slurm scripts are self-contained -- arguments are baked into the
slurm file. Variant changes are made by editing the slurm script, not
via `--export`.

## Run order

Run `build_library.slurm` first, then `subsample_baseline.slurm` (it
reads the library characteristics produced by the first step). Stage 04
`baseline_comparison.py` consumes both sets of outputs to render Fig 7.

## Outputs

```
outputs/03_kirsch_library/
  build_library/<slug>/{config.json, library.npy,
                        characteristics.json, characteristics.npz}
  subsample_baseline/{lhs,sobol,random}/<slug>/
                       {config.json, subsample_<method>.json}
```

Slug for `build_library` encodes `n{n_traces}_t{n_years}_ssi{ssi}_s{seed}`;
slug for `subsample_baseline` encodes the source library slug plus
`_n{n_select}_s{seed}`.

## Figures

None at the stage level. Fig 7 is rendered by
[workflows/04_moea_find_single_site/plots/baseline_comparison.py](../04_moea_find_single_site/plots/baseline_comparison.py)
and finalized for the manuscript by
[workflows/99_manuscript_figures/make_figures.py](../99_manuscript_figures/make_figures.py).
