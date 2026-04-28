# Stage 08 — NYC reservoir sensitivity

## Purpose

Apply common global sensitivity analysis (SA) methods to the structured
drought-hazard realization sample produced by MOEA-FIND, quantifying how
NYC reservoir outcomes (combined-storage drawdown, minimum storage
fraction, Montague flow reliability and vulnerability) respond to each
optimized drought-hazard characteristic. SD (Stage 07) and SA (Stage 08)
are complementary applications of well-established exploratory-modeling
methods to drought-hazard space; the contribution is the *application
to hazard space*, not the methods themselves. Backs SI-12.

The factor space is exactly the set of MOEA-FIND objective axes that
drove the upstream realization sample, read from the upstream
`results.json::objective_keys`. Non-optimized drought characteristics
are excluded by design — only the optimized dimensions inherit the
structured (epsilon-tiled, near-uniform) coverage that justifies SA in
the first place.

Three sample-free SA methods are run on every outcome:

- **Delta moment-independent (Borgonovo, 2007)** — manuscript anchor.
- **PAWN density-based (Pianosi & Wagener, 2015)** — comparator.
- **RBD-FAST (Tarantola et al., 2006)** — comparator.

Sobol / Morris / FAST are excluded: they require dedicated Saltelli /
radial / frequency-coded sampling designs and would force a re-run of
Pywr-DRB on a different sample, which defeats the point of re-using
the existing Stage-06 metric bank.

## Run order

Must run **after** stage 06 (consumes
`outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet` and
the upstream MOEA-FIND `results.json`). All downstream of this stage is
figure production. Stage 07 runs in parallel — both consume the same
metric bank but neither depends on the other.

## Drivers

| Driver | SLURM | YAML config | Outputs | Working figures | Section |
|---|---|---|---|---|---|
| [run_sa.py](run_sa.py)                 | [slurm/run_sa.slurm](slurm/run_sa.slurm) | [configs/all_methods.yaml](configs/all_methods.yaml) | `outputs/exp10_nyc_sensitivity/<slug>/`         | `outputs/exp10_*/<slug>/figures/` | §SI-12 |
| [compare_methods.py](compare_methods.py) | (no SLURM; serial)                    | (consumes Stage-08 outputs)                          | `figures/nyc_sensitivity/comparisons/`          | `figures/nyc_sensitivity/comparisons/` | §SI-12 |

`run_sa.py` accepts `--config <path>` to load a YAML preset; CLI flags
still override. Four presets in [configs/](configs/) cover the headline
all-methods run plus single-method re-runs.

`compare_methods.py` reads multiple `run_sa.py` output dirs and emits
side-by-side tornadoes plus run × run rank-correlation heatmaps. Use it
for sensitivity-of-method studies (different upstream metric sets,
different sample sizes, different SA hyperparameters).

## Inputs

- **X (factors).** Per-realization drought characteristics from the
  upstream MOEA-FIND archive (`pareto_chars` field in
  `outputs/exp04_kirsch_single_site/<slug>/results.json`).
- **Y (outcomes).** Per-realization NYC outcomes from the Stage-06
  metric bank (`outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet`).

Realization ids must align between the two files. The driver inner-joins
on the index and drops rows with non-finite values.

## Outputs

Per Stage-08 run under `outputs/exp10_nyc_sensitivity/<slug>/`:

```
config.json                                # invocation record
results/
  indices_<method>.parquet                 # long-form: outcome × factor × indices
  bootstrap_<method>.parquet               # bootstrap CI per (outcome, factor)
  rank_stability_<method>.parquet          # bootstrap-rank Spearman summary
  convergence_<method>.parquet             # index vs n
  cross_method_rank_corr.parquet           # method × method (per outcome)
  cross_outcome_rank_corr.parquet          # outcome × outcome (per method)
  selection_log.json                       # criterion result + chosen anchor
figures/                                   # only with --plot
  tornado_<outcome>_<method>.pdf
  heatmap_indices_<method>.pdf             # factor × outcome
  convergence_<outcome>.pdf                # method-faceted
  cross_method_rank_corr_<outcome>.pdf
  cross_outcome_rank_corr_<method>.pdf
```

The slug encodes the upstream archive, the metric set, the methods, the
outcome count, and the seed — see [src/slugs.py](../../src/slugs.py).

## Manuscript-method selection criterion

Computed per outcome, applied post-HPC. For an SA method to qualify as
the manuscript anchor it must **simultaneously**:

1. Bootstrap CI on the top-ranked factor not span zero.
2. Bootstrap-rank Spearman ≥ `rank_spearman_threshold` (default 0.8).
3. Cross-method rank correlation ≥ `cross_method_threshold` (default
   0.7) with at least one other method on the same outcome.

If Delta passes, it anchors the figure for that outcome. If Delta fails,
PAWN is checked next, then RBD-FAST. The decisions and the failing
condition (if any) are recorded in `selection_log.json`. The criterion
is fixed before HPC to avoid post-hoc method shopping.

## Manuscript figures

- Main: 1–2 figures (post-HPC choice — candidate A is a tornado for the
  selected outcome under the selected method; candidate B is a factor ×
  outcome heatmap under the selected method).
- SI: §SI-12 — method comparison (PAWN vs Delta vs RBD-FAST tornadoes
  per outcome), cross-method rank-correlation matrices, cross-outcome
  rank-correlation matrices, sample-size convergence curves.

Working PDFs land in `outputs/exp10_*/<slug>/figures/`. Final assets
are regenerated by `99_manuscript_figures/make_figures.py`.

## SI figure budget

§SI-12 targets 2–3 multi-panel figures: method comparison + cross-method
rank correlation, cross-outcome rank correlation, sample-size
convergence. Sized to fit alongside SI-11 (scenario-discovery
manifold).

## Dependencies

- `SALib >= 1.5.2` (numpy 2.0 compatible) is a core dependency declared
  in [pyproject.toml](../../pyproject.toml). It is pulled in automatically
  by `pip install -e .`.
- `scipy.stats.spearmanr` (already a transitive dependency).
- Standard runtime: numpy, pandas, matplotlib, pyarrow.
