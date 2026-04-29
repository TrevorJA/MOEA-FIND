# MOEA-FIND Publication Plan

*Living document. Phase tracker, exit criteria, and research questions for the
single-site Cannonsville methods paper targeted at *Water Resources Research*.*

The authoritative section outline lives in `drafts/manuscript_main_draft.md`.
The current snapshot of source code and workflows lives in `code_state.md`.
Section-level prose decisions live in `governance/design_decisions.md` (DD-01
through DD-N). Reviewer-anticipation work lives in
`reference/reviewer_defenses.md`.

Scope is a methods paper with a single-site Cannonsville case study. The
multi-site Delaware River Basin extension is deferred to a follow-up paper and
is not a main-text result. Terminology is "drought hazard space"; banned
substitutes are listed in `governance/style_guide.md` §5.11–5.17.

The Phase labels in this document (A–E) are workflow phases for pre-HPC, HPC,
and submission blocks. They are not the HPC compute-phase labels (beta, gamma,
delta) used elsewhere.

---

## 1. Research Questions

The paper is organised around three primary questions and three supporting
questions. Items resolved or superseded since the original RQ list have been
folded into design decisions and are noted with their DD pointer.

**RQ1. Structured coverage via multi-objective search.** Can a multi-objective
evolutionary algorithm with an L1 Manhattan-norm auxiliary objective produce
synthetic streamflow ensembles with structured coverage of a user-defined
drought characteristic space? Evaluation: L2-star discrepancy and
nearest-neighbour CV against LHS, Sobol, and library-LHS baselines (Section
3.2). Resolved on the analytic benchmark for K=2..6 (DD-11); pending on
Cannonsville (DD-12 audit).

**RQ2. Dimensionality scaling and interior coverage.** Does the construction
deliver interior-filling coverage of the K-dimensional feasible image as K
grows? **Resolved (DD-11) on the analytic K-ball benchmark at K=2..6.** The
high-dimensional Cannonsville analogue is open (DD-12). See
`evidence/shell_vs_interior_diagnostic.md`.

**RQ3. Physical plausibility.** Do generated traces satisfy the plausibility
envelope of the Kirsch-Nowak generator (autocorrelation, seasonal cycle,
flow-duration curve, cross-site correlation if multi-site)? Constraints are
calibrated per-experiment by bootstrap (`workflows/0N_<stage>/diag_constraint_calibration.py`)
and reported per-run in SI; see DD-05 and DD-14 for the production constraint
choice.

**RQ4. Comparison to alternatives.** How does MOEA-FIND compare to library
LHS-subsampling and to repeated single-objective FIND runs in coverage and in
NFE budget? Evaluation in Section 3.2 of the manuscript with the Cannonsville
single-site result.

**RQ5. Operational scenario discovery.** Does structured hazard-space coverage
yield a better classifier of an *operational* failure label (Pywr-DRB storage
or flow-target breach) than equal-size LHS-in-decision-space ensembles? See
Part A of `experiment_plan.md` for the redesigned Section 3.3 protocol; see
`reference/reviewer_defenses.md` HC-2 for the prior circular-label issue this
redesign addresses.

**RQ6. Epsilon and anti-ideal sensitivity.** How do epsilon values and
anti-ideal placement affect ensemble size, spacing, and feasible-region
coverage? See DD-04 for anti-ideal placement and `code_alignment_backlog.md`
Item 8 for tightening D* in the production Cannonsville run.

Folded or deferred:

- *RQ on event-level vs trace-level framing* — covered by DD-01.
- *RQ on multi-site basin application* — deferred to follow-up paper; not a
  main-text result for this submission.
- *RQ on BART scenario discovery* — replaced by gradient boosted tree
  classification in Section 3.3 per the experiment plan; BART is no longer
  in scope.

---

## 2. Phase Tracker

Each task lists exit criteria. Do not advance to the next phase until prior
exit criteria are satisfied or have been explicitly deferred in writing.

### Phase B. Experiment hardening (mostly local; partially complete)

| ID | Task | Exit criterion |
|----|------|----------------|
| B1 | Epsilon × NFE sensitivity sweep on the analytic benchmark, deferred to HPC. Driver: `workflows/01_analytic_validation/eps_nfe_sweep.py`. | Coverage-vs-epsilon curve produced; recommended default ε stored alongside the aggregate output. |
| B2 | Single-site Kirsch + SSI-3 objectives + Manhattan, NFE ≥ 50 000. Driver: `workflows/04_moea_find_single_site/run_moea_find.py` in `residual` mode. | Pareto front exists; plausibility spot-check passes (lag-1 ACF, FDC, seasonal cycle). |
| B3 | Re-run B2 with bootstrap-calibrated constraints. | Constrained Pareto retains ≥80 % of unconstrained hypervolume; no plausibility failures. |
| B4 | Generate 10 000-trace Kirsch library and subsample (LHS, Sobol, random). Drivers: `workflows/0N_<stage>/05_kirsch_library_build.py`, `workflows/0N_<stage>/06_library_subsample_baseline.py`. | Comparison table and headline coverage figure produced. |
| B5 | Plausibility report (ACF, FDC, seasonal cycle, Hurst). Currently bundled into script 04 via `--plot`. | Diagnostic figure in `figures/`. |

### Phase C. HPC transition

| ID | Task | Exit criterion |
|----|------|----------------|
| C1 | MM Borg Python wrapper installed on HPC; smoke test on a standard MOEA test problem. | Wrapper runs end-to-end on login or interactive node. |
| C2 | Wire MM Borg into `workflows/04_moea_find_single_site/run_moea_find.py` with checkpointing (`passNFE_ALH_PyCheckpoint` branch). | Checkpoint-restart test passes (identical Pareto within ε after resume). |
| C3 | Cannonsville single-site MOEA-FIND runs at production NFE and ≥5 seeds. | Seed-averaged Pareto archive; coverage metrics reported vs library subsample. |
| C4 | Section 3.3 Pywr-DRB policy re-evaluation on the Cannonsville archive and on the LHS baseline. Driver: `workflows/0N_<stage>/09_drb_policy_reeval.py`. | Operational failure labels computed; classifier comparison produced. |

### Phase D. Manuscript drafting (parallel with Phase C)

| ID | Task | Exit criterion |
|----|------|----------------|
| D1 | Sections 1 (Introduction) and 2 (Methods) prose. | ≥1500 words; every cited paper verified against Zotero; style guide passes. |
| D2 | Section 3 (Results) paired with the seven main-text figures. | Each figure has a paragraph stating the quantitative result before interpretation. |
| D3 | Sections 4 (Discussion) and 5 (Conclusions). | Limitations explicitly cite DD-10 (feasible-region framing) and the Kirsch historical-envelope constraint. |
| D4 | Supporting Information populated. | Every SI item cross-referenced from main text. |
| D5 | Internal review pass; co-author review. | Tracked changes incorporated; open comments resolved or recorded. |

### Phase E. Reproducibility and submission

| ID | Task | Exit criterion |
|----|------|----------------|
| E1 | Each experiment reproducible from a single `sbatch` (HPC) or `python` (local smoke test) invocation. Pinned seeds in the SLURM wrappers. | Fresh clone reproduces all figures with one command per experiment. |
| E2 | `REPRODUCE.md` written. | Third-party reviewer can rerun core results. |
| E3 | DOI-mint code release on Zenodo via GitHub release. | DOI cited in manuscript. |
| E4 | Submit preprint (EarthArXiv) simultaneously with journal submission. | Submission receipts archived. |

---

## 3. Checkpoint gates

**Gate 1 — End of Phase B.**
Single-site Kirsch MOEA-FIND runs reproducibly to ≥50 000 NFE. Plausibility
report passes. Library baseline (≥10 k traces) archived. Coverage comparison
figure tells a coherent story whichever direction the result goes.

**Gate 2 — End of Phase C.**
MM Borg runs on HPC with checkpointing verified. Cannonsville single-site
ensemble generated with ≥5 seeds. All main-text figures generated at draft
quality. Section 3.3 operational labels computed.

**Gate 3 — End of Phase D.**
Internal review pass complete. Every claim has a script, figure, or citation
backing it. Limitations section acknowledges DD-01, DD-10, DD-12, and the
Kirsch historical-envelope constraint.

---

## 4. Deferred / out of scope

- Comparison to NSGA-III, MOEA/D, or other MOEAs on either the analytic
  benchmark or the production Cannonsville problem. MM Borg MOEA is the
  only optimizer used in the project per DD-07 (updated 2026-04-28); the
  earlier EpsNSGAII (platypus) stand-in for laptop runs has been removed.
- Phase randomization and KDE-smoothed CDF generators (DD-02 options E, F):
  mentioned in Discussion only.
- Multi-site Delaware River Basin application: deferred to a follow-up paper.
- Demand as a scenario-discovery dimension: deferred per
  `code_alignment_backlog.md` Item 9.
- Generalisation to k=4,5 drought objectives in the applied study: limited to
  the analytic K-ball scaling figure.
