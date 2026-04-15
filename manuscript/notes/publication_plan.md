# MOEA-FIND Publication Plan and Task List

*Created: 2026-04-13 (Session 15). Living document.*

> **2026-04-14 sync update.** This file predates the finalised
> manuscript structure. The authoritative current state of the
> manuscript is `manuscript_main_draft.md` and the active plan is
> `C:\Users\tjame\.claude\plans\staged-painting-dongarra.md`. Key
> alignment points: (a) the manuscript scope is a methods paper
> with a single-site Cannonsville case study, not a
> Delaware-River-Basin case-study paper; the multi-site extension
> is now explicitly deferred to a follow-up paper and is not a
> main-text result; (b) the manuscript uses a five-section flat
> outline (Introduction, Methods, Results, Discussion, Conclusions)
> with no Introduction subheadings and seven main-text figures;
> (c) the terminology is "drought hazard space"; "outcome space",
> "feasibility discovery", "admissible drought characteristics",
> and "Reed group" are disallowed in body text per style guide
> sections 5.11 through 5.17; (d) the scenario discovery
> demonstration in Section 3.3 uses gradient boosted trees, not
> BART or PRIM; (e) the shell-versus-interior question has been
> resolved empirically across $K = 2$ to $K = 6$ as DD-11, and
> DD-12 opens the high-dimensional Cannonsville empirical audit.
> Phase labels in this document (Phase A through Phase E) remain
> valid as workflow phases for the pre-HPC, HPC, and submission
> blocks and are not the same as the HPC phase labels (beta,
> gamma, delta) used elsewhere.

This document consolidates the state of the MOEA-FIND project and defines the work required to bring it to publication-level completeness. It is scoped to a single methods paper (**target venue: *Water Resources Research***) with the Cannonsville single-site study as the applied case and the multi-site Delaware River Basin demonstration deferred to follow-up. Phases are ordered so that each checkpoint produces a self-contained, reviewable artifact.

---

## 1. Current State Assessment (as of 2026-04-13)

### 1.1 Completed

- **Conceptual foundation.** Manhattan-norm hyperplane trick documented and justified (`notes/method_proposal.md`, §2). Design decisions DD-01 through DD-10 resolved or deferred with rationale. Literature review reconciled against Zotero; fabricated precedent claims removed (DD-09).
- **Core library.** `src/objectives.py` (drought metrics, SSI, Manhattan norm), `src/constraints.py`, `src/analysis.py` (discrepancy, NN metrics), `src/kirsch_wrapper.py` (SynHydro integration with index + residual modes), `src/library.py` (library-and-subsample baseline), `src/parametric.py` (kappa4 + D-vine), `src/data.py`, and `src/experiment_utils.py` are implemented. Unit tests cover objectives, constraints, analysis, and the Kirsch wrapper.
- **Phase 1 validation.** Analytic 2D and 3D experiments complete. Hyperplane constraint satisfied to machine precision (std ~1e-16). Coverage metrics computed against LHS / Sobol / Random baselines. See DD-10 for interpretation.
- **Figures 1-6.** Publication-quality PNG+PDF for the conceptual and analytic-validation figures are generated.
- **Baseline scaffolding.** `scripts/05_kirsch_library_build.py` and `scripts/06_library_subsample_baseline.py` implement the library generation and subsampling comparison targets; ready to run at scale on HPC.
- **Repository reorganization (2026-04-13).** All numerical experiments now live under `scripts/` as a flat, numbered layout (`scripts/NN_<slug>.py` + `scripts/NN_<slug>.slurm`). Each script is pinned to a specific manuscript section and figure via `scripts/README.md`. `_common.sh` centralizes SLURM module loads, venv activation, and MPI launch selection. The previous `experiments/` tree has been replaced by a redirect README.
- **Plotting consolidation (2026-04-13).** All plotting functions live under `src/plotting/` (`analytic.py`, `coverage.py`, `drought_space.py`, `trace_diagnostics.py`, `convergence.py`, `style.py`). Each function is tagged with the manuscript section and figure it produces. Experiment scripts call these via `--plot`; `scripts/10_plot_manuscript_figures.py` is the single-command assembler that reads every `outputs/expNN_*/` and (re)generates the complete publication figure set in `figures/`. Deleted plot scripts from the old `experiments/` layout were harvested for reusable panels (Figs 1-4 from `plot_poc_results.py`, correlation and spacing panels from `plot_diagnostics.py`) and integrated into the `src/plotting/` modules.
- **Outputs tree (2026-04-13).** `outputs/` now follows a strict `outputs/expNN_<slug>/` convention with `config.json` + `results.json` + optional `pareto.npz` + optional `cells/` for SLURM arrays. Legacy results from the pre-reorganization layout have been moved under `outputs/archive/` and are read-only. `.gitignore` keeps `outputs/README.md` tracked while ignoring all regenerable artifacts. See `outputs/README.md` for the full convention.

### 1.2 In progress / blocked

- **Phase 2 Kirsch experiments.** Wrapper and scripts exist; no archived single-site run has been executed past smoke-test NFE budgets. Blocked only on compute time and a chosen epsilon grid.
- **Parametric experiments (Phase 3).** Module exists; no run has been validated end-to-end.
- **DRB multi-site application (Phase 4).** Only a placeholder README.
- **HPC transition.** MM Borg Python wrapper and checkpointing resources identified (Feb/Aug 2025 Water Programming posts, `passNFE_ALH_PyCheckpoint` branch) but not yet wired up.
- **Manuscript.** `manuscript/manuscript_main_draft.md` and `manuscript/supporting_info_draft.md` are empty stubs. This plan bootstraps their section structure.

### 1.3 Known gaps and risks

1. **Apples-to-apples coverage comparison.** DD-09 Option D (library subsampling in drought space) is the only fair baseline; needs a 10 000+ trace Kirsch library generated before the paper's headline comparison figure can be produced.
2. **Feasible-region framing.** DD-10 shows that raw uniformity is *not* MOEA-FIND's headline advantage. The paper must reframe around *structured coverage of the physically feasible drought region* or readers will compare NN_CV to Sobol and conclude the method underperforms.
3. **Event-level vs. trace-level framing (DD-01).** The paper promises both. Only the trace-level framing is implemented in code.
4. **Plausibility audit.** No end-to-end trace plausibility report exists (autocorrelation, FDC, seasonal cycle, Hurst). Required before reviewers will accept the bootstrap pipeline.
5. **Borg vs. platypus.** Current experiments use serial Borg via the SynHydro-adjacent pipeline; no confirmation that the Borg binary path / license is set up on HPC.

---

## 2. Work Breakdown by Phase

Each task lists **exit criteria** that define "done" for the purpose of a checkpoint review. Do not advance to the next phase until the prior phase's exit criteria are satisfied (or have been explicitly deferred in writing).

### Phase A. Repository and manuscript scaffolding (local, ~1 week)

| ID | Task | Exit criterion |
|----|------|----------------|
| A1 | Populate `manuscript/manuscript_main_draft.md` with section headings, one-sentence abstracts per section, and placeholder figure callouts. | A reader can infer the paper's narrative arc from the outline alone. |
| A2 | Populate `manuscript/supporting_info_draft.md` mirroring the SI sections in `notes/research_questions.md`. | SI sections named, each with scope note. |
| A3 | Create `notes/publication_plan.md` (this file) and link from `README.md` and `CLAUDE.md`. | Plan discoverable from project root. |
| A4 | Freeze a `v0.1-poc` git tag capturing the state of Phase 1 results. | Tag exists locally. Push deferred to user. |
| A5 | Audit `.gitignore` for data, binaries, Borg sources, `outputs/`, `figures/*.png`. Confirm no secrets or licensed Borg code has ever been committed. | `git log -- lib/borg*` empty; `outputs/` ignored. |

### Phase B. Experiment hardening (local, ~2 weeks)

| ID | Task | Exit criterion |
|----|------|----------------|
| B1 | **Exp 1.3 / script 03.** Epsilon × NFE sensitivity sweep on 3D analytic problem. **Deferred to HPC (2026-04-13):** local attempt stopped after the first few cells because platypus archive sorting + the O(n²) L2-star discrepancy calculation together make the 20k–50k NFE cells impractical on a workstation. Driver `scripts/03_eps_nfe_sweep.py` is HPC-ready (`--mode cell`, `--task-id` for SLURM array, `--metric-cap 0` to disable subsampling). Launcher `scripts/03_eps_nfe_sweep.slurm` is a 90-task SLURM array. Local partial log archived at `outputs/exp_1_3/sweep_log_partial_local.txt`. | Coverage-vs-epsilon curve produced. Recommended default ε stored in `outputs/exp03_eps_nfe_sweep/aggregate.json` and documented in `notes/method_proposal.md` §4. |
| B2 | **Exp 2.1 / script 04.** Single-site Kirsch + SSI-3 objectives (duration, avg severity) + Manhattan. NFE ≥ 50 000. `scripts/04_kirsch_single_site.py` in `residual` mode. | Pareto front exists; plausibility spot-check passes (lag-1 acf within 0.05 of historical; no zero-flow months; FDC envelope overlap > 90 %). |
| B3 | **Exp 2.2 / script 04 `--constrained`.** Re-run with lag-1 autocorrelation + non-drought annual statistics constraints. Report constraint-violation rate at each NFE decile. | Constrained Pareto has ≥80 % of unconstrained hypervolume and no plausibility failures. |
| B4 | **Exp 2.3 / scripts 05 + 06.** Generate 10 000-trace Kirsch library via `scripts/05_kirsch_library_build.py` (MPI) and subsample via `scripts/06_library_subsample_baseline.py` (LHS, Sobol, random). Headline Figure 7 comparison. | Comparison table + figure produced. Narrative updated to match whichever result actually holds. |
| B5 | **Exp 2.4 / script 07.** Event-level short-trace formulation (5–10 yr DVs, event-level metrics: duration, peak, cumulative severity, onset month). Requires event-level helpers in `src/objectives.py` which are not yet implemented; `scripts/07_event_level_kirsch.py` currently runs only under `--dry-run`. | At least one event-level Pareto archive and plausibility check. |
| B6 | **Plausibility report.** Script producing Fig 6 analogue: acf, FDC, seasonal cycle, Hurst for bootstrap and (future) parametric traces. Currently bundled into script 04 via `--plot`; may split into `06a_plausibility_report.py` if it grows. | PDF committed to `figures/`. |

### Phase C. HPC transition (≈1 week of wall time)

| ID | Task | Exit criterion |
|----|------|----------------|
| C0 | **Exp 1.3 HPC sweep.** `sbatch scripts/03_eps_nfe_sweep.slurm` launches the 90-task SLURM array (6 ε × 3 NFE × 5 seeds). Follow with `python scripts/03_eps_nfe_sweep.py --mode aggregate`. Doubles as HPC env smoke test. | `outputs/exp03_eps_nfe_sweep/aggregate.json` written; recommendation stored in `notes/method_proposal.md` §4. |
| C1 | Build MM Borg Python wrapper per Feb 2025 Water Programming post on target HPC. Smoke test with DTLZ2. Update `CLUSTER_*` variables in `scripts/_common.sh` to match the cluster's module and partition names. | Wrapper runs end-to-end on login/interactive node. |
| C2 | Port `KirschBorgWrapper` problem definition to MM Borg with checkpointing (`passNFE_ALH_PyCheckpoint` branch). Wire into `scripts/04_kirsch_single_site.py` under an MPI code path. | Checkpoint-restart test passes (identical Pareto within epsilon after resume). |
| C3 | **Script 08.** Run Exp 4.1 (multi-site, 4 DRB inflow sites, shared indices). `sbatch scripts/08_drb_multisite_moea.slurm` (8-seed SLURM array). NFE ≥ 500 000 per seed. | Seed-averaged Pareto archive; cross-site correlation matrix reported. |
| C4 | **Script 08** (continued). Run Exp 4.2 with 3 drought objectives + Manhattan. | Structured ensemble saved. Coverage metrics reported vs. library subsample. |
| C5 | **Script 09.** Stage Exp 4.3 policy re-evaluation hand-off to `NYCOptimization`. | Ensemble exported in the format Pywr-DRB expects; `scripts/09_drb_policy_reeval.py` wired up. |

### Phase D. Manuscript drafting (parallel with Phase C, ~3 weeks)

| ID | Task | Exit criterion |
|----|------|----------------|
| D1 | Draft Introduction + Background/Related Work. | ≥ 1500 words; every cited paper verified against Zotero. |
| D2 | Draft Methods (formulation, wrapper, objectives, constraints, comparison baselines). | Reproduces notation in the figures. |
| D3 | Draft Results sections paired with Figs 1-7. | Each figure has a paragraph that states the quantitative result before interpretation. |
| D4 | Draft DRB case-study section with Figs 8-9. | Gated on Phase C completion. |
| D5 | Draft Discussion, Limitations, Conclusions. | Limitations explicitly cite DD-10 feasible-region framing and Kirsch historical-envelope constraint. |
| D6 | Populate SI §§ SI-1 through SI-7. | Every SI item cross-referenced from the main text. |
| D7 | Internal review pass (code + manuscript), then co-author review. | Tracked changes incorporated; open comments resolved or recorded. |

### Phase E. Reproducibility and submission (≈1 week)

| ID | Task | Exit criterion |
|----|------|----------------|
| E1 | Make every experiment reproducible from a single `sbatch scripts/NN_*.slurm` invocation (HPC) or `python scripts/NN_*.py` invocation (local smoke test). Pinned seeds live in the `.slurm` wrappers; any per-run overrides are captured in the per-run `outputs/expNN_*/config.json`. | Fresh clone reproduces all figures with one command per experiment. |
| E2 | Write `REPRODUCE.md` describing HPC and local paths, data acquisition, Borg license handling. | Third-party reviewer can rerun core results. |
| E3 | DOI-mint the code release (Zenodo via GitHub release). | DOI cited in the manuscript. |
| E4 | Submit preprint (EarthArXiv) simultaneously with journal submission. | Submission receipts archived in `notes/archive/`. |

---

## 3. Checkpoint Criteria (Go / No-Go Gates)

Use these gates before declaring a phase "done." Each gate is intentionally conservative; failing a gate means looping back, not papering over.

**Gate 1 — End of Phase B (ready for HPC).**
- Single-site Kirsch MOEA-FIND runs reproducibly to ≥50 000 NFE.
- Plausibility report shows autocorrelation, FDC, and seasonal cycle errors within published tolerances for stochastic streamflow generators.
- Library baseline (≥10k traces) exists and is archived.
- Exp 2.3 comparison figure tells a coherent story, *whichever direction* the result goes.
- Manuscript outline (A1, A2) in place.

**Gate 2 — End of Phase C (ready for results writing).**
- MM Borg runs on HPC with checkpointing verified.
- DRB multi-site ensemble generated with ≥8 independent seeds.
- Cross-site correlation errors < 0.05 relative to historical.
- All figures referenced in §5 of `research_questions.md` are generated at draft quality.

**Gate 3 — End of Phase D (ready for submission).**
- Internal review pass complete.
- Every claim in the manuscript has a script, figure, or citation backing it.
- Limitations section acknowledges DD-10, DD-01, and the Kirsch historical-envelope constraint.

---

## 4. Recommended Next Steps (immediately actionable)

1. **Finalize manuscript scaffolding (A1, A2) today.** Empty files are the biggest friction against coherent writing sessions. The outlines in `manuscript/manuscript_main_draft.md` and `manuscript/supporting_info_draft.md` should be the first thing collaborators see.
2. **Commit the publication plan + manuscript scaffolds, then tag `v0.1-poc`.** Gives a stable reference point to diff against during Phase B.
3. **HPC access.** Local attempt at Exp 1.3 (2026-04-13) confirmed that the ε×NFE sweep is impractical on a workstation. File the HPC access request immediately; Exp 1.3 is now the first job queued on it (Phase C0), and it doubles as the smoke test for the HPC Python + platypus/Borg environment.
4. **Start the 10 000-trace Kirsch library build** (Exp 2.3) on HPC as soon as access lands. This is the bottleneck for the headline coverage figure.
5. **Begin single-site Exp 2.1 locally** at ≤ 20k NFE with the Kirsch wrapper so that the plausibility tooling (acf, FDC, Hurst) is validated end-to-end before HPC scale-up.
6. **Reframe the narrative.** Before writing, decide explicitly: is the paper's headline "uniform coverage" (which DD-10 weakens) or "structured coverage of the feasible drought region" (which DD-10 supports)? The Phase B results will confirm, but the manuscript outline should commit to the latter now.

---

## 5. Deferred / out of scope

- Comparison to NSGA-III / MOEA/D / platypus epsilon-MOEA (user directive: Borg only).
- Phase randomization and KDE-smoothed CDF tracks (DD-02 options E, F) — mentioned in discussion only.
- Integration with BART scenario discovery beyond Exp 4.3 hand-off; that is a sequel paper.
- Generalization to k=4,5 drought objectives in the applied study; limited to scaling figure on the analytic problem.
