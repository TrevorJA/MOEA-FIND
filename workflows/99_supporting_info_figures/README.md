# 99_supporting_info_figures

Every **supporting-information** figure script. Each script calls the
reusable plotting functions in `src/plotting/` and reads cached upstream
`outputs/<stage>/.../` artifacts — it never re-runs an experiment.

Created by the spring-cleaning reorg: the former
`workflows/<stage>/plots/*.py` SI/diagnostic drivers were consolidated
here so there are no per-stage plotting folders or SLURM jobs.

## Run

One job builds them all (best-effort; a script whose upstream is absent
is logged and skipped, not fatal):

```
sbatch workflows/99_supporting_info_figures/slurm/make_si_figures.slurm [SLUG]
```

`SLUG` is passed as `--slug` to scripts that accept it. Per-script slug
overrides can be refined inside the slurm as the SI matures.

## Contents

~36 scripts spanning SI-1 (analytic hyperplane), SI-A..D calibration
diagnostics, SI-F/G/H stage-04 ablations & bounded-metric mapping
checks, SI-I coverage verification, and SI-12/12b NYC sensitivity.
Main-text figures are built separately by
`workflows/99_manuscript_figures/`.
