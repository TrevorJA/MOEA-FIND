# Stage 01 -- Analytic validation


## Purpose

Demonstrate the L1 + epsilon-tile MOEA-FIND device on a synthetic
objective space whose true Pareto manifold is known analytically. Backs
the manuscript's analytic-validation section and the K=2..6 dimension
sweep that motivates applying MOEA-FIND to constrained hydrologic
problems where uniform space filling is non-trivial.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [analytic_2d.py](analytic_2d.py)         | K=2 visual proof: L1 + tile yields uniform Pareto.        | Fig 1 |
| [analytic_3d.py](analytic_3d.py)         | K=3 with hypercube feasible region.                       | Fig 2 |
| [eps_nfe_sweep.py](eps_nfe_sweep.py)     | Joint sensitivity of coverage to (epsilon, NFE).          | Fig 3 |
| [dimension_sweep.py](dimension_sweep.py) | Coverage of K=1..6 hypercube under fixed epsilon.         | Fig 4 |

`dimension_sweep.py` defaults to the **hypercube** feasibility shape; the
K-ball variant remains available via `--feasible-shape ball` for
historical comparison but is no longer the manuscript default.

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/01_analytic_validation/<driver>/<slug>/`. Figures are produced
by the paired plotting drivers under
[plots/](plots/) which read those artifacts and write to
`figures/01_analytic_validation/<driver>/<slug>/`. Re-rendering a
figure never requires re-running an experiment.

## Slurm

| Compute slurm | Plotting slurm |
|---|---|

`eps_nfe_sweep` is a SLURM array (4 epsilons x 3 NFEs x 3 seeds = 36
cells, one core each, fits a single 40-core node). After the array
completes, run the aggregator:

```bash
python workflows/01_analytic_validation/eps_nfe_sweep.py --mode aggregate
```

then submit the plotting slurm.

`dimension_sweep` is a 5-task array over K=2..6 (one core each).

All slurm scripts are self-contained -- arguments are baked into the
slurm file. Variant changes are made by editing the slurm script, not
via `--export`.

## Outputs

```
outputs/01_analytic_validation/
  analytic_2d/<slug>/{config.json, results.json, pareto.npz}
  analytic_3d/<slug>/{config.json, results.json, pareto.npz}
  eps_nfe_sweep/{config.json, cells/cell_<id>.json, aggregate.json}
  dimension_sweep/k{2..6}/{config.json, results.json, samples.npz}
```

## Figures

```
figures/01_analytic_validation/
  analytic_2d/<slug>/analytic_2d.pdf
  analytic_3d/<slug>/analytic_3d.pdf
  eps_nfe_sweep/eps_nfe_sweep.pdf
  dimension_sweep/dimension_sweep.pdf
```

Final manuscript PDFs are regenerated from these working figures by
[workflows/99_manuscript_figures/make_figures.py](../99_manuscript_figures/make_figures.py).
