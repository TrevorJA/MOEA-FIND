# Stage 06 — Pywr-DRB policy re-evaluation

## Purpose

Run Pywr-DRB simulation against each Pareto drought produced by stages
04 or 05, producing a per-realization metric bank (FFMP exposure, NYC
storage, flow targets, delivery reliability). This stage is **simulation
only** — post-processing (satisficing labelling, scenario discovery)
lives in stage 07.

## Run order

Must run **after** stage 04 (single-site) or stage 05 (multi-site)
because it consumes a Pareto archive. `verify_drought_coverage.py` runs
**before** the simulation as a sanity check on the input Pareto's
drought-space coverage. The simulation output `metric_bank.parquet`
feeds stage 07.

## Drivers

| Driver | SLURM | Outputs | Working figures | SI subsection |
|---|---|---|---|---|
| [verify_drought_coverage.py](verify_drought_coverage.py) | [slurm/verify_drought_coverage.slurm](slurm/verify_drought_coverage.slurm) | `outputs/diag_verify_drought_coverage/<slug>/` | `outputs/diag_*/<slug>/figures/` | §SI-I |
| [policy_reeval.py](policy_reeval.py)                     | [slurm/policy_reeval.slurm](slurm/policy_reeval.slurm)                     | `outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet` | `outputs/exp09_*/<slug>/figures/` | §7.3 main |

`policy_reeval.py` uses `src.pywrdrb_bridge` to (1) replay Pareto DVs
to multi-site daily, (2) preprocess via pywrdrb's
`PredictedInflowEnsemblePreprocessor`, (3) batch-simulate across
realizations, and (4) compute the metric bank via
`src.satisficing_metrics.compute_metric_bank`. Expensive — can take
hours on full Pareto archives.

## Manuscript figures

- Main: contributes data to Fig 9 (satisficing map; produced by stage
  07's `scenario_discovery_plots.py`).
- SI: §SI-I (drought-coverage verification of the Pareto inputs).

## SI figure budget

§SI-I targets 2-3 multi-panel figures verifying that the Pareto archive
covers the drought-characteristic space without gaps before Pywr-DRB
runs are committed.
