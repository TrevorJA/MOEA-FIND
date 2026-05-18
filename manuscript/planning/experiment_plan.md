# Experiment Plan — Sections 3.2, 3.3 and 3.4

*Originally merged 2026-04-15 from `section3_3_redesign.md` and
`hpc_deployment_status.md`. Trimmed 2026-04-27 to align with the current
`workflows/` structure and to retire HC-2-resolved content. The current
authoritative code state lives in `code_state.md`.*

---

## Part A. Section 3.3 — operational failure labels (HC-2 resolved)

The original Section 3.3 design defined a binary failure label as the joint
condition that mean drought severity $D_1$ and mean drought duration $D_2$
both exceed their respective historical 80th percentiles, then trained a
gradient boosted tree classifier on $(D_1, D_2)$ as features. The
circularity is exact: the failure threshold is defined on the same
coordinates that MOEA-FIND optimises over, so any space-filling sample in
$(D_1, D_2)$ produces a better-calibrated classifier of a threshold in
$(D_1, D_2)$ than a density-weighted library subsample does. This was
reviewer critique HC-2 in `reference/reviewer_defenses.md` and is now
resolved by replacing the threshold label with a Pywr-DRB-derived
operational outcome.

**Production protocol.** Two equal-size ensembles enter Pywr-DRB simulation:
the MOEA-FIND archive and an equal-size library LHS subsample drawn in
generator decision-variable space (Bonham et al., 2024). Each member
provides Cannonsville inflows; other sites draw from the historical
Kirsch-Nowak baseline. Operating rules follow the standard FFMP
implementation in the current Pywr-DRB model. Demand is held at the
average historical NYC Delaware Basin per-capita value. Initial conditions
use the median observed Cannonsville storage on 1 April. Simulation length
is matched to trace length.

**Failure label.** Three Hashimoto et al. (1982) metrics are computed.
Supply reliability $\lambda$ is the fraction of months in which combined
NYC storage stays above the FFMP Level 1 trigger. Vulnerability $\nu$ is
the maximum single-month shortfall as a fraction of demand. Flow target
reliability $\mu$ is the fraction of months in which Montague and Trenton
non-NYC flow targets are both met. The binary failure label $y = 1$ if
$\lambda < 0.95$ or $\nu > 0.15$, with thresholds calibrated against the
historical record and threshold sensitivity reported in Supporting
Information.

**Classifier.** A gradient boosted tree classifier (Friedman, 2001;
Chen and Guestrin, 2016) is trained on each ensemble separately to predict
$y$ from $(D_1, D_2, D_3)$. Hyperparameters are tuned by five-fold
cross-validation and held fixed across both ensembles. Comparison metrics
are AUC on a held-out test set, Brier score, and decision boundary
visualisation.

The full pre-resolution discussion lives in `archive/critical_analysis_v1.md`
and in `reference/reviewer_defenses.md` HC-2.

---

## Part B. Production pipeline (current code)

Part B in earlier versions of this document tracked TODOs T-01 through T-10
against an `scripts/hpc/section32_*.slurm` layout that has since been
superseded. The current production pipeline lives under
`workflows/0N_<stage>/` with one-to-one SLURM wrappers under
`workflows/0N_<stage>/slurm/`. The relevant scripts are:

| Stage | Script | Purpose |
|-------|--------|---------|
| §3.2 single-site MOEA | `workflows/04_moea_find_single_site/run_moea_find.py` | Single-site Kirsch + SSI-3 objectives + Manhattan; production-ready. |
| §3.2 library | `workflows/0N_<stage>/05_kirsch_library_build.py` | 10 k+ Kirsch-Nowak library generation. |
| §3.2 baseline | `workflows/0N_<stage>/06_library_subsample_baseline.py` | LHS / Sobol / random subsampling of the library. |
| §3.2 figures | `workflows/0N_<stage>/10_plot_manuscript_figures.py` | Headline coverage figure assembly. |
| §3.2 baseline comparison | `workflows/0N_<stage>/11_baseline_comparison.py` | MOEA-FIND vs library coverage. |
| §3.3 policy re-eval | `workflows/0N_<stage>/09_drb_policy_reeval.py` | Four-stage Pywr-DRB pipeline: replay Pareto DVs → multi-site daily → simulate → classify. Uses `src/pywrdrb_bridge.py`, `src/satisficing_metrics.py`, `src/satisficing_labels.py`, `src/scenario_discovery.py`. |
| §3.3 satisficing sweep | `workflows/0N_<stage>/12_satisficing_sweep.py` | Sweep satisficing thresholds for robustness. |
| §3.3 plots | `workflows/0N_<stage>/17_drb_scenario_discovery_plots.py` | Decision boundary and feature-importance figures. |

Diagnostics that gate the pipeline:

| Diagnostic | Purpose |
|-----------|---------|
| `workflows/0N_<stage>/diag_constraint_calibration.py` | Bootstrap-calibrate the five hydrologic constraints per trace length. |
| `workflows/0N_<stage>/diag_dv_uniformity_calibration.py` | Calibrate the DV-uniformity ablation constraint. |
| `workflows/0N_<stage>/diag_kirsch_convergence.py` | Pareto size and hypervolume vs NFE. |
| `workflows/0N_<stage>/diag_kirsch_wrapper_fidelity.py` | Validate Kirsch wrapper output against direct SynHydro. |
| `workflows/0N_<stage>/diag_kirsch_wrapper_geometry.py` | Geometric properties of the wrapper DV mapping. |
| `workflows/0N_<stage>/diag_shell_vs_interior.py` | Shell-vs-interior K=2..6 sweep underlying DD-11 and `evidence/shell_vs_interior_diagnostic.md`. |

For the current snapshot of `src/` and `workflows/` and the alignment between
each module and the manuscript, see `code_state.md`. For per-experiment
TODOs, decisions, and open methodology choices, see `code_alignment_backlog.md`
and `governance/design_decisions.md`.

---

## Part C. DD-15 — joint metric-and-T justification protocol

*Added 2026-04-29 as the gating prerequisite for Part B production
runs. Until DD-15 produces (K*, T*), all §3.2 single-site MOEA-FIND
results in this document and in the draft manuscript are PROVISIONAL.*

DD-15 (`governance/design_decisions.md`) defines a five-stage SLURM
pipeline that resolves DD-01 (trace length T) and DD-04 (drought-metric
set K) jointly, replacing the current "by analogy to prior art"
defenses with empirical evidence.

**T grid (coarse-bracket):** `{5, 10, 20, 30}` water years.
**K cardinality:** {3, 4} jointly evaluated.
**Reference data:** USGS 01423000 Cannonsville inflow, 73 water years
1950-10-01 to 2023-09-30, single-site (matches the §3.2 production
case study).
**Output of record:** `outputs/02_calibration/decision_matrix/pareto_front_KxT.json` —
recommended (K*, T*) plus two alternates plus the cost-vs-T fit.

| Stage | Driver | Outputs |
|-------|--------|---------|
| 1 — historical T-blocks | `workflows/02_calibration/t_sensitivity_historical.py` (array, one task per T) + `t_sensitivity_aggregate.py` | per-T spread, degeneracy, correlations, K-set alternatives; cross-T summary, cluster-stability ARI, K-set Jaccard |
| 2 — Kirsch fidelity | `workflows/03_kirsch_library/build_library_extended.py` (array, one task per T) + `t_sensitivity_kirsch_compare.py` | 10 k baseline Kirsch realizations × 4 Ts × 28 metrics; KS, Frobenius corr-shift, L2-star coverage |
| 3 — joint K × T decision | `workflows/02_calibration/eval_cost_timing.py` (array) + `decision_matrix.py` | per-T MOEA inner-loop wall-time; 5-component composite score; (K*, T*) recommendation |
| 4 — confirmatory MM Borg | `workflows/04_moea_find_single_site/run_moea_find.py` with `--n-years T*` and `--metric-set <K*-preset>` | Pareto archive at (K*, T*); 120 ranks × 8 h |
| 5 — manuscript figures | `src/plotting/02_calibration/figure_{a,b,c,d}_*.py` | Figure A (metric stability across T), B (Spearman clustermap preservation), C (drought-space coverage), D (exemplar traces) |

The pipeline is chained with `--dependency=afterok:` from
`workflows/02_calibration/run_t_sensitivity.sh` (orchestrator
blueprint). User explicitly approves the (K*, T*) recommendation
before Stage 4 commits 120-rank Borg compute (≈ 960 core-hours).
Stage-1 viability gate revised 2026-04-29 to drop the saturation-as-
degeneracy criterion: saturation is informational only (extreme-event
metrics like `max_duration` saturate at long T but remain handled by
the spread-screen IQR>0 check, which is independent of T).

**Carry-through invariants.** All stages preserve DD-11 (Manhattan
formulation: `f_j = D_j + ‖D − D*‖_1`, never `|D_j − D*_j|`), DD-14
(AD DV-uniformity constraint), and the SSI-3/SSI-12/Q80 lock-in (all
calibrators fitted on the full historical record once and reused
across T). Multi-site (DD-02) remains deferred.

---

## Part D. Section 3.4 — magnitude-varying sensitivity in drought-hazard space (DD-16)

The diagnostic complement to §3.3 scenario discovery. §3.3 asks
whether operational failure is *separable* in the
drought-characteristic space; §3.4 asks which characteristic
dimensions *dominate* at each severity of an operational hazard
outcome, adapting Hadjimichael et al. (2020) magnitude-varying
sensitivity analysis (MV-SA) with the factor space and magnitude
axis both inhabiting the empirical drought-hazard / outcome space
that MOEA-FIND structures coverage over.

**Code-side status (2026-04-30).** Engine, plot library, CLI driver,
SLURM scripts, and synthetic-data tests are landed:
`src/magnitude_varying_sa.py`,
`src/plotting/magnitude_varying_sa.py`,
`workflows/09_magnitude_varying_sa/`,
`tests/test_magnitude_varying_sa.py`. All 8 sanity tests pass on
synthetic data (severity-dependent ranking recovered; control factor
sits below real factors at every percentile; degenerate slices
return NaN without crash).

**Production runs.**

1. **Stage A — proof-of-concept on existing K=3 archive.** ✓ COMPLETE
   (job 217983, 2026-04-30, 8 min wall). Response form:
   `within_trace_percentile` (Hadjimichael Variant 2); Delta method
   only; N=3308 full archive; n_bootstrap=50; 19 ranks (1 per
   percentile), single-node exclusive. Output slug:
   `mvsa__src=residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad__axis=nyc_min_storage_frac__resp=within_trace_percentile__metric_set=h30e96e__methods=delta__n_perc=19__s=42`.
   Key result: rank crossover at τ ≈ 0.40; `mean_magnitude` leads at
   severe percentiles, `time_in_drought_fraction` leads at
   moderate/favorable; control Δ ≈ 0.05–0.09 throughout.

2. **Stage B — figures.** ✓ COMPLETE (job 217984, 2026-04-30). PDFs at
   `figures/09_magnitude_varying_sa/run_mv_sa/<slug>/`:
   `stacked_area_delta.pdf`, `lines_with_ci_delta.pdf`.

3. **Stage C — re-run on K* archive.** Once DD-15 Stage 4
   confirmatory MOEA-FIND run lands, re-run Stages A and B against
   that slug with n_bootstrap=200 and all three methods (Delta, PAWN,
   RBD-FAST). One-line invocation; no methodology change.

4. **Stage D — SI-12b robustness panels.** Repeat Stages A–B with
   `--magnitude-axis montague_flow_vulnerability` and
   `--magnitude-axis nyc_drawdown_days_below_0.25` for SI-12b
   robustness. Same archive, different invocations.

5. **Stage E — conditional-response variant.** Repeat Stage A with
   `--response-form conditional --secondary days_L6` (FFMP Level 6
   trigger from the Pywr-DRB drought operating rule) for a
   supporting SI-12b diagnostic.

**Gating.** Stages A–B complete. Stages C–E unblocked except Stage C
which is gated on DD-15 completion. The methodology in
`src/magnitude_varying_sa.py` is K-agnostic.

**Carry-through invariants.** Same factor-space rule as Stage 08
(read `objective_keys` from the upstream archive's
`results.json`; non-optimised characteristics excluded). Bootstrap
seed 42; percentile grid frozen at 19 evenly-spaced points; control
uniform factor included in every run. MPI: 19 ranks (one per
percentile), single-node exclusive, NO mpi4py collective calls —
each rank loads data independently and writes a partial parquet +
sentinel file; rank 0 polls sentinels and concatenates.
