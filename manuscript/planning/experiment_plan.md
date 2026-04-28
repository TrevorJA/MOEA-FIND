# Experiment Plan — Sections 3.2 and 3.3

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
