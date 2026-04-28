# Stage 08 — NYC reservoir sensitivity (placeholder)

## Status

**Not yet implemented.** This stage is a placeholder for the planned
exploratory-modeling sensitivity analysis of NYC reservoir outcomes
under the MOEA-FIND drought ensembles. The implementation is deferred
to a later session and will be specified during a separate planning
discussion.

## Purpose (planned)

Quantify how NYC reservoir outcomes (storage depletion, FFMP exposure,
delivery reliability, deficit volume) depend on drought-hazard
characteristics (duration, severity, frequency, peak intensity).
Treat the MOEA-FIND Pareto front as a structured sample over the
drought hazard space and use exploratory-modeling tools (PRIM, GBT,
random forest sensitivity, partial-dependence) to map outcome
distributions to characteristic axes.

## Run order (planned)

Will consume:

- Pareto droughts from stage 04 (single-site) or stage 05 (multi-site)
- Pywr-DRB outputs from stage 06 (`metric_bank.parquet`)
- Optionally satisficing labels from stage 07

Will not write back into upstream stages.

## Drivers (planned)

To be defined.

## Manuscript figures (planned)

- Main: TBD (likely 1-2 sensitivity-driven scenario-discovery figures).
- SI: TBD.

Working figures will land in `figures/nyc_sensitivity/`.
