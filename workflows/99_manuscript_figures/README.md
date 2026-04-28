# Stage 99 — Manuscript figure regeneration

## Purpose

Single source of truth for the final manuscript figure set. Reads every
upstream `outputs/expNN_*/` directory and regenerates each main-text
and SI figure into `figures/main/` and `figures/supplementary/`. The
*only* script allowed to write into those two folders.

## Run order

Runs **last**, after every contributing upstream stage has been
executed. Subsets can be regenerated via `--only fig0X` for fast
iteration.

## Drivers

| Driver | SLURM | Outputs |
|---|---|---|
| [make_figures.py](make_figures.py) | [slurm/make_figures.slurm](slurm/make_figures.slurm) | `figures/main/fig*.pdf`, `figures/supplementary/figSI*.pdf` |

`make_figures.py` does not run any expensive computation — it is a
plotting-only aggregator that imports figure factories from
`src.plotting.*` and saves the resulting PDFs. If a contributing
upstream output is missing, the corresponding figure is skipped with a
warning rather than failing the run.

## Manuscript figures

Every figure that appears in the manuscript main text or supplementary
material is produced here. The promotion contract:

- Main-text figures: `figures/main/figNN_<descriptor>.pdf`
- SI figures: `figures/supplementary/figSI_<descriptor>.pdf`
- Working/exploratory figures from individual stages: themed subfolders
  (`figures/analytic/`, `figures/calibration/`, ...). Promotion to
  `main/` or `supplementary/` happens here.

## CLI

```bash
python workflows/99_manuscript_figures/make_figures.py            # regenerate all
python workflows/99_manuscript_figures/make_figures.py --only fig04  # one figure
```
