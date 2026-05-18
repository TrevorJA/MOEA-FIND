# Supporting Information: MOEA-FIND

*SCAFFOLD DRAFT. Section headers only. Each subsection body is a placeholder
with a pointer to the underlying authoritative source. Prose will be
drafted only after the corresponding main-text section has been written
against final empirical results, so that the SI is written to support
specific main-text claims rather than anticipating them.*

*Companion to the main text. Each section is cross-referenced from a
specific location in the main text. Target length is approximately three
thousand words plus figures and tables.*

*Last updated: 2026-04-28 (workflow reorganization; stage→subsection mapping added).*

## Workflow stage to SI subsection mapping

Each method-consideration driver under [workflows/](../../workflows/) earns
SI text. The mapping below is the contract: every kept driver in stages 02
and 04, plus the dimension sweep in stage 01 and the verification/discovery
drivers in stages 06–07, is referenced below. Drivers that produce only
intermediate calibration data (no SI figures) are still tracked, with the
subsection that consumes them noted.

| Driver | SI subsection |
|---|---|
| `01_analytic_validation/eps_nfe_sweep.py`            | SI-3 (epsilon and NFE sensitivity) |
| `01_analytic_validation/dimension_sweep.py`          | SI-2 (interior-filling coverage) |
| `02_calibration/constraint_calibration.py`           | SI-4b (constraint regime ablation, hydrologic arm) |
| `02_calibration/dv_uniformity_calibration.py`        | SI-4b (constraint regime ablation, DV-space arm) |
| `02_calibration/wrapper_fidelity.py`                 | SI-5 (DV parameterisation, fidelity check) |
| `02_calibration/wrapper_geometry.py`                 | SI-5 (DV parameterisation, geometry check) |
| `02_calibration/kirsch_convergence.py`               | SI-4 (Borg convergence + Kirsch wall-clock) |
| `02_calibration/metric_blocks.py`                    | SI-6 (legacy single-T metric blocks; superseded by DD-15) |
| `02_calibration/metric_explorer.py`                  | SI-6 (legacy single-T screening; uses `src/metric_screening.py` shared API) |
| `02_calibration/t_sensitivity_historical.py`         | SI-6 (Stage-1 per-T historical-block screening, DD-15) |
| `02_calibration/t_sensitivity_aggregate.py`          | SI-6 (Stage-1 cross-T aggregation, DD-15) |
| `03_kirsch_library/build_library_extended.py`        | SI-6 (Stage-2 28-metric Kirsch ensemble, DD-15) |
| `02_calibration/t_sensitivity_kirsch_compare.py`     | SI-6 (Stage-2 KS / Frobenius / coverage, DD-15) |
| `02_calibration/eval_cost_timing.py`                 | SI-6 (Stage-3 per-T eval cost regression, DD-15) |
| `02_calibration/decision_matrix.py`                  | SI-6 (Stage-3 joint K × T matrix + (K*, T*), DD-15) |
| `04_moea_find_single_site/wrapper_mode_ablation.py`  | SI-5 (DV parameterisation ablation) |
| `04_moea_find_single_site/wrapper_mode_compare.py`   | SI-5 (DV parameterisation ablation) |
| `04_moea_find_single_site/dv_uniformity_ablation.py` | SI-4b (constraint regime ablation) |
| `04_moea_find_single_site/dv_uniformity_compare.py`  | SI-4b (constraint regime ablation) |
| `04_moea_find_single_site/event_level.py`            | SI-9 (event-level Kirsch — placeholder below) |
| `06_pywrdrb_reeval/verify_drought_coverage.py`       | SI-10 (Pareto coverage verification — placeholder below) |
| `07_scenario_discovery/satisficing_sweep.py`         | SI-11 (satisficing manifold + GBT — placeholder below) |
| `07_scenario_discovery/scenario_discovery_plots.py`  | SI-11 (satisficing manifold + GBT) |
| `08_nyc_sensitivity/run_sa.py`                       | SI-12a (NYC aggregate sensitivity to drought-hazard characteristics) |
| `08_nyc_sensitivity/compare_methods.py`              | SI-12a (NYC aggregate sensitivity, cross-run comparison) |
| `09_magnitude_varying_sa/run_mv_sa.py`               | SI-12b (NYC magnitude-varying sensitivity to drought-hazard characteristics; main-text §3.4 anchor) |

Each subsection targets approximately 2-3 multi-panel figures (clean,
multi-panel academic style) plus supporting tables. Subsections SI-9,
SI-10, and SI-11 below are new placeholders added with this mapping.

---

## SI-1. Manhattan-distance auxiliary objective and the K-dimensional coverage argument

*Cross-referenced from main-text Section 2.2.3.*

> *Placeholder. The L1 formulation is fixed: $f_j = D_j$ for $j = 1, \ldots,
> K$ and $f_{K+1} = \lVert D - D^* \rVert_1$. The non-dominance argument and
> the codimension-one affine-subspace identity are documented in
> `governance/design_decisions.md` §DD-11 and implemented in
> `src/objectives.py`. The formal derivation will be drafted into this
> section once the empirical sections that motivate the framing exist.*

---

## SI-2. Empirical verification of interior-filling coverage

*Cross-referenced from main-text Section 3.1.*

> *Placeholder. The dimension-sweep evidence on a constrained analytic
> benchmark is documented in `evidence/shell_vs_interior_diagnostic.md`
> (MM Borg MOEA on HPC, per DD-07). The benchmark, dimensionalities
> tested, and final tabulation will be finalised once §2.3 of the main
> text is written.*

---

## SI-3. Epsilon and function-evaluation sensitivity on the analytic test problem

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. Driver: `workflows/01_analytic_validation/eps_nfe_sweep.py`. SLURM
> array configuration in `workflows/0N_<stage>/slurm/`. Recommended default epsilon
> vector pending HPC sweep completion.*

---

## SI-4. Convergence diagnostics for Borg MOEA

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. Multi-seed hypervolume curves, archive size traces, and
> adaptive operator selection probabilities for every main-text experiment.
> Pending HPC production runs.*

---

## SI-4b. Constraint regime ablation

*Cross-referenced from main-text Section 2.2.4 and `governance/design_decisions.md` §DD-14.*

> *Placeholder. Constraint regime ablation. The set of regimes compared and
> the production choice are not yet final. Current drivers and diagnostics
> live in `workflows/0N_<stage>/` and `workflows/0N_<stage>/`; see
> `planning/code_state.md`.*

---

## SI-5. Sensitivity to the decision-variable parameterisation

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. The decision-variable parameterisation of the Borg–generator
> coupling is not yet final. The production choice and any wrapper-mode
> ablation will be drafted here once settled. Current drivers live in
> `workflows/0N_<stage>/`; see `planning/code_state.md`.*

---

## SI-6. Joint sensitivity to drought metric set and trace length T (DD-15)

*Cross-referenced from main-text Sections 2.1 and 2.4.*

> *Placeholder for the four DD-15 figures plus the decision-matrix
> table. Once Stage 3 produces
> `outputs/02_calibration/decision_matrix/pareto_front_KxT.json`, this
> subsection will report:*
>
> 1. **Stage-1 distributional stability across T.** Per-metric
>    descriptors (median, IQR, robust spread, skew) on stride-1
>    historical T-blocks at T ∈ {5, 10, 20, 30}. Adjusted Rand Index
>    of cluster memberships across T-pairs and Jaccard score of
>    top-K sets vs the T=20 reference. **Key small-T finding:** the
>    SSI-12 family of metrics zero-degenerates in ~43% of T=5 blocks
>    (no SSI-12 events fit in 60 months) but the Tier-A/C/D/E/F
>    candidate pool still produces hundreds of feasible strict-rung
>    K=3 sets. **Key long-T finding:** at T=30, every block contains
>    the historical-worst drought event, so `max_duration` and
>    `worst_severity` saturate at their record values (informational;
>    handled by the IQR>0 spread screen, which excludes them
>    automatically when they cease to vary).
> 2. **Stage-2 Kirsch fidelity at every T.** Per-metric Kolmogorov–
>    Smirnov statistic comparing the 10 000-realization baseline-
>    Kirsch distribution to the historical T-block distribution.
>    Frobenius norm of the elementwise Spearman-correlation
>    difference. L2-star discrepancy and NN-CV in each candidate
>    K-set's subspace.
> 3. **Stage-3 K × T decision matrix.** Five-component composite
>    score (min robust spread × constraint-rung penalty × 1−max|ρ|
>    × 1−KS_max × 1−norm cost) for the top-5 K=3 and top-5 K=4
>    strict-rung sets at every surviving T. Per-T MOEA inner-loop
>    median wall-time and projected 120-rank wall-time at NFE=200 000
>    (preliminary timings as of 2026-04-29: T=5 → 92 ms/eval, T=10 →
>    145 ms, T=20 → 248 ms, T=30 → 343 ms).
> 4. **Figures A–D.** Figure A — metric stability across T (per-metric
>    violin + Kirsch ridge + MOEA scatter, columns = T, rows = K* +
>    runners-up). Figure B — three-panel Spearman clustermap (Hist,
>    Kirsch, MOEA) at T*. Figure C — pairwise scatter matrix in K*
>    space with L2-star + NN-CV side bar. Figure D — three exemplar
>    MOEA-FIND traces with SSI-3 shading and Hist-IQR annotation.
>
> *Driver pointers:*
> [`workflows/02_calibration/t_sensitivity_historical.py`](../../workflows/02_calibration/t_sensitivity_historical.py),
> [`t_sensitivity_aggregate.py`](../../workflows/02_calibration/t_sensitivity_aggregate.py),
> [`build_library_extended.py`](../../workflows/03_kirsch_library/build_library_extended.py),
> [`t_sensitivity_kirsch_compare.py`](../../workflows/02_calibration/t_sensitivity_kirsch_compare.py),
> [`eval_cost_timing.py`](../../workflows/02_calibration/eval_cost_timing.py),
> [`decision_matrix.py`](../../workflows/02_calibration/decision_matrix.py),
> and the four Stage-5 figure scripts under
> [`src/plotting/02_calibration/`](../../src/plotting/02_calibration/).
> Reference: `governance/design_decisions.md` §DD-15;
> `planning/experiment_plan.md` Part C.

---

## SI-7. Per-site plausibility diagnostics (multi-site extension)

*Cross-referenced from main-text Section 2.4. Included only if the
multi-site extension makes the main text.*

> *Placeholder. Multi-site DRB extension is currently deferred to a
> follow-up paper per `planning/publication_plan.md` §4. This section will
> most likely be cut from the single-site submission. Site list and
> diagnostics are not finalised.*

---

## SI-8. Pareto archive reference tables

*Cross-referenced from main-text Section 3.2.*

> *Placeholder. Reference tables for the production archive will be
> released as supplementary files at submission. The exact columns released
> will be defined once §3.2 is written.*

---

## SI-9. Event-level Kirsch objective formulation

*Cross-referenced from main-text Section 2.2 (Fig 5 inset).*

> *Placeholder. Event-level objective formulation that ties Kirsch DVs
> directly to drought-event characteristics (rather than aggregate
> series-level statistics). Driver:
> [`workflows/04_moea_find_single_site/event_level.py`](../../workflows/04_moea_find_single_site/event_level.py).
> Target 2-3 multi-panel figures: per-event Pareto, sensitivity to event
> definition window, comparison to series-level objectives. Pending HPC
> production runs.*

---

## SI-10. Pareto-archive drought coverage verification

*Cross-referenced from main-text Section 3 (Pywr-DRB re-evaluation
inputs).*

> *Placeholder. Verifies that the Pareto archive supplied to Pywr-DRB
> covers the drought-characteristic space without gaps before the
> expensive simulation is committed. Driver:
> [`workflows/06_pywrdrb_reeval/verify_drought_coverage.py`](../../workflows/06_pywrdrb_reeval/verify_drought_coverage.py).
> Target 2-3 multi-panel figures: drought-space coverage scatter,
> per-axis marginal histograms, FDC subset diagnostics.*

---

## SI-11. Satisficing manifold and GBT classifier diagnostics

*Cross-referenced from main-text Section 3.3 / Fig 9.*

> *Placeholder. Per-definition GBT decision-boundary overlays in
> drought-feature space, classifier ROC summary, manifest status table.
> Drivers:
> [`workflows/07_scenario_discovery/satisficing_sweep.py`](../../workflows/07_scenario_discovery/satisficing_sweep.py)
> and
> [`workflows/07_scenario_discovery/scenario_discovery_plots.py`](../../workflows/07_scenario_discovery/scenario_discovery_plots.py).
> Manifest of binary satisficing rules:
> [`workflows/07_scenario_discovery/satisficing_manifest.yaml`](../../workflows/07_scenario_discovery/satisficing_manifest.yaml).
> Target 2-3 multi-panel figures.*

---

## SI-12. NYC sensitivity to drought-hazard characteristics

*Cross-referenced from main-text Sections 3.3 (scenario discovery)
and 3.4 (magnitude-varying sensitivity).*

### SI-12a. Aggregate sensitivity (Stage 08)

> *Placeholder. Applies common global sensitivity analysis methods
> (Delta moment-independent as the manuscript anchor; PAWN density-based
> and RBD-FAST as comparators) to the structured drought-hazard
> realization sample produced by MOEA-FIND. Factor space is exactly the
> set of MOEA-FIND objective axes that drove the upstream archive
> (read from `results.json::objective_keys`); non-optimized drought
> characteristics are excluded by design because they have no space-
> filling guarantee. Outcomes carried through SA: NYC minimum combined
> storage fraction, NYC drawdown days below 25%, Montague flow
> reliability, Montague flow vulnerability — main-text outcome chosen
> post-HPC. Driver:
> [`workflows/08_nyc_sensitivity/run_sa.py`](../../workflows/08_nyc_sensitivity/run_sa.py).
> Cross-run comparison driver:
> [`workflows/08_nyc_sensitivity/compare_methods.py`](../../workflows/08_nyc_sensitivity/compare_methods.py).
> Target figure budget: 2-3 multi-panel SI figures (method comparison
> tornadoes, cross-method and cross-outcome rank-correlation matrices,
> sample-size convergence). Pending HPC production runs.*

### SI-12b. Magnitude-varying sensitivity (Stage 09)

> *Placeholder. Adapts Hadjimichael et al. (2020) magnitude-varying
> sensitivity analysis to the drought-hazard factor space (the §3.4
> headline). For each percentile τ ∈ {0.05, 0.10, …, 0.95} of an
> operational hazard outcome, sensitivity indices for every
> drought-hazard characteristic (plus a uniform-random control
> factor) are computed and presented as a stacked-area panel.
> Magnitude axes carried through MV-SA: NYC minimum combined storage
> fraction (main-text §3.4 anchor), NYC drawdown days below 25%, and
> Montague flow vulnerability (SI robustness panels). Conditional-
> response variant on FFMP Level 6 (Pywr-DRB drought operating-rule
> trigger) is reported as a supporting diagnostic. Triangulation
> across Delta moment-independent, PAWN density-based, and RBD-FAST
> establishes method consensus; bootstrap CIs (n=200, seed-pinned)
> bound each percentile slice. The control factor establishes the
> magnitude-varying noise floor. Driver:
> [`workflows/09_magnitude_varying_sa/run_mv_sa.py`](../../workflows/09_magnitude_varying_sa/run_mv_sa.py);
> plotting driver:
> [`src/plotting/09_magnitude_varying_sa/run_mv_sa.py`](../../src/plotting/09_magnitude_varying_sa/run_mv_sa.py).
> Target figure budget: 3 multi-panel SI figures (per-axis stacked
> area with method panels; per-axis line plots with bootstrap
> ribbons; conditional-response variant). Pending HPC production
> run; proof-of-concept uses the existing K=3 archive
> `residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad`. Re-run
> on the post-DD-15 K* archive once Stage 4 lands.*

---

> *Figure and table numbering in this Supporting Information will be
> finalised once the main-text figures are locked. Any SI section whose
> upstream output has not been produced by the submission date will be cut
> or replaced with a data-tables-only entry.*
