# Stage 08 -- NYC sensitivity analysis


## Purpose

Apply common global sensitivity analysis (SA) methods to the structured
drought-hazard realization sample produced by MOEA-FIND, quantifying how
NYC reservoir outcomes (combined-storage drawdown, minimum storage
fraction, Montague flow reliability and vulnerability) respond to each
optimized drought-hazard characteristic. Backs SI-12.

The factor space is exactly the set of MOEA-FIND objective axes that
drove the upstream realization sample, read from the upstream
`results.json::objective_keys`. Non-optimized drought characteristics
are excluded by design -- only the optimized dimensions inherit the
structured (epsilon-tiled, near-uniform) coverage that justifies SA in
the first place.

Three sample-free SA methods are run on every outcome:

- **Delta moment-independent (Borgonovo, 2007)** -- candidate anchor.
- **PAWN density-based (Pianosi & Wagener, 2015)** -- comparator.
- **RBD-FAST (Tarantola et al., 2006)** -- comparator.

Sobol / Morris / FAST are excluded: they require dedicated Saltelli /
radial / frequency-coded sampling designs and would force a re-run of
Pywr-DRB on a different sample, which defeats the point of re-using
the existing Stage-06 metric bank.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [run_sa.py](run_sa.py)                      | Compute SA indices, bootstrap CIs, rank-stability, cross-method/outcome rank correlations, and the per-outcome anchor-selection log. | SI-12 |
| [plots/run_sa.py](plots/run_sa.py)          | Render tornado, heatmap, convergence, and rank-correlation diagnostics from cached parquets.                                          | SI-12 |
| [plots/compare_methods.py](plots/compare_methods.py) | Side-by-side comparison of factor rankings across `run_sa` slugs (different upstream metric sets / sample sizes / hyperparameters). | SI-12 |

`run_sa.py` accepts `--config <path>` to load a YAML preset; CLI flags
still override. Presets in [configs/](configs/) cover the headline
all-methods run plus single-method re-runs.

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/08_nyc_sensitivity/<driver>/<slug>/`. Figures are produced by
the centralized figure jobs in workflows/99_*_figures/ which read those
artifacts and write to `figures/08_nyc_sensitivity/<driver>/<slug>/`.
Re-rendering a figure never requires re-running SA.

## Run order

Must run **after** Stage 06 (consumes
`outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet`
and the upstream MOEA-FIND
`outputs/04_moea_find_single_site/run_moea_find/<src_slug>/results.json`).
All downstream of this stage is figure production. Stage 07 runs in
parallel -- both consume the same metric bank but neither depends on
the other.

## Slurm

| Compute slurm | Plotting slurm |
|---|---|

All slurm scripts are self-contained -- arguments are baked into the
slurm file. Variant changes are made by editing the slurm script, not
via `--export`.

## Inputs

- **X (factors).** Per-realization drought characteristics from the
  upstream MOEA-FIND archive (`pareto_chars` field in
  `outputs/04_moea_find_single_site/run_moea_find/<src_slug>/results.json`).
- **Y (outcomes).** Per-realization NYC outcomes from the Stage-06
  metric bank
  (`outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet`).

Realization ids must align between the two files. The driver inner-joins
on the index and drops rows with non-finite values.

## Outputs

```
outputs/08_nyc_sensitivity/
  run_sa/<slug>/
    config.json                              # invocation record
    results/
      indices_<method>.parquet               # long-form: outcome x factor x indices
      bootstrap_<method>.parquet             # bootstrap CI per (outcome, factor)
      rank_stability_<method>.parquet        # bootstrap-rank Spearman per (outcome, factor)
      rank_stability_summary_<method>.parquet
      convergence_<method>.parquet           # index vs n
      cross_method_rank_corr.parquet         # method x method (per outcome)
      cross_outcome_rank_corr.parquet        # outcome x outcome (per method)
      selection_log.json                     # criterion result + chosen anchor
```

The slug encodes the upstream archive, the metric set, the methods, the
outcome count, and the seed -- see [src/slugs.py](../../src/slugs.py).

## Figures

```
figures/08_nyc_sensitivity/
  run_sa/<slug>/
    tornado_<outcome>_<method>.pdf
    heatmap_indices_<method>.pdf             # factor x outcome
    convergence_<outcome>.pdf                # method-faceted
    cross_method_rank_corr_<outcome>.pdf
    cross_outcome_rank_corr_<method>.pdf
  compare_methods/<tag>/
    compare_tornado_<method>_<outcome>.pdf
    compare_rho_<method>_<outcome>.pdf
```

Final manuscript PDFs are regenerated from these working figures by
[workflows/99_manuscript_figures/make_figures.py](../99_manuscript_figures/make_figures.py).

## Manuscript-method selection criterion

Computed per outcome, applied post-HPC (decisions persisted to
`selection_log.json`). For an SA method to qualify as the manuscript
anchor it must **simultaneously**:

1. Bootstrap CI on the top-ranked factor not span zero.
2. Bootstrap-rank Spearman >= `rank_spearman_threshold` (default 0.8).
3. Cross-method rank correlation >= `cross_method_threshold` (default
   0.7) with at least one other method on the same outcome.

If Delta passes, it anchors the figure for that outcome. If Delta fails,
PAWN is checked next, then RBD-FAST. The decisions and the failing
condition (if any) are recorded in `selection_log.json`. The criterion
is fixed before HPC to avoid post-hoc method shopping.

## Dependencies

- `SALib >= 1.5.2` (numpy 2.0 compatible) is a core dependency declared
  in [pyproject.toml](../../pyproject.toml). It is pulled in
  automatically by `pip install -e .`.
- `scipy.stats.spearmanr` (already a transitive dependency).
- Standard runtime: numpy, pandas, matplotlib, pyarrow.
