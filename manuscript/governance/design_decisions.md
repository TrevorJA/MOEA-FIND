# Design Decisions

*Living document. One entry per decision. Older versions of any entry can be
recovered from git history; resolved debates and pre-decision option tables
are not preserved here. When a decision is reversed or superseded, mark the
old entry SUPERSEDED with a pointer and write the new entry as a separate DD.*

**Status conventions.** Every entry carries one of:

- **SETTLED** — implemented and not expected to change.
- **PENDING** — investigation underway; resolution expected.
- **OPEN** — question not yet investigated; may change anything downstream.
- **SUPERSEDED** — preserved for provenance only.

The active gating item is DD-12 (Cannonsville empirical audit; PENDING).
DD-04 (drought metrics) was settled 2026-04-27 with the ``primary``
preset (see entry below).

---

## DD-01: Trace length and event-vs-trace framing

*Status: PENDING — superseding numerical defense underway under DD-15
(2026-04-29). The 2026-04-12 entry below is preserved for context until
DD-15 produces (T*, K*).*

The production case study is single-site Cannonsville with moderate-length
traces (on the order of decades), matching the trace-level framing of
Borgomeo et al. (2015), Zaniolo et al. (2024), and Wheeler et al. (2025).
Aggregate drought metrics computed across the full trace.

The historical default `n_years_out = 20` in
`src/experiment_config.py:ExperimentConfig.n_years_out` was inherited by
analogy to that prior art and never numerically defended. Short-trace
framings (T ∈ {5, 10}) are also under consideration on the grounds that
short blocks better correspond to *distinct hazard periods*, strengthening
the hazard-space interpretation of MOEA outcomes and distinguishing this
work from the "few long series" prevailing standard.

An event-level framing in which short traces (1–8 years) carry a single
designed drought event remains a follow-up direction. Event-level
metrics (peak intensity, cumulative severity, recovery rate) are natural
substitutes for trace-level aggregates and would shrink the decision-variable
count substantially. Not the production focus.

The exact trace length and the trace-vs-event question are tied to the
choice of drought metrics in DD-04 and resolve jointly under DD-15.

**Pointer:** `drafts/manuscript_main_draft.md` §2.4;
`workflows/04_moea_find_single_site/run_moea_find.py`; DD-15.

---

## DD-02: Decision-variable formulation

*Status: SETTLED for production (2026-04-12).*

Production runs use a continuous decision-variable vector in $[0, 1]^d$,
converted by the wrapper into either historical year indices or empirical-CDF
residuals (the two injection modes implemented in `src/kirsch_wrapper.py`;
see DD-06). The continuous formulation is compatible with Borg's variational
operators and is the only mode used by the production driver.

A parametric track (Kappa-4 marginals plus D-vine copula for temporal
dependence) is a POC in `src/parametric.py` and is not integrated into the
production pipeline. It is a follow-up direction motivated by the bootstrap
generator's inability to extrapolate beyond the historical envelope.

**Pointer:** `src/kirsch_wrapper.py`; `src/parametric.py` (POC); DD-06.

---

## DD-03: Bootstrap granularity

*Status: SETTLED 2026-04-13.*

The Kirsch pipeline uses $B = 1$ (independent monthly resampling) followed
by Cholesky correlation imposition in normal-score space, with a shifted-year
construction for the December–January boundary. Mechanical block bootstraps
($B > 1$) are not used; their within-block-correlation benefit is replaced
by the Cholesky step, and they introduce block-boundary discontinuities.

**Pointer:** SynHydro Kirsch generator; DD-06.

---

## DD-04: Objective set (drought metrics)

*Status: PROVISIONAL 2026-04-27 (production default for development);
joint numerical justification with trace length T underway under DD-15
(2026-04-29). Final K* tuple may differ from the `primary` preset
described below pending Stage 3 decision matrix.*

The production drought metric set is the ``primary`` preset of
:mod:`src.drought_metrics`: three continuous, literature-anchored
metrics that span depth, volume, and persistence axes of drought:

1. **Mean event severity** (``mean_severity``) — mean of the per-event
   minimum SSI value (depth axis). Continuous; standard event-level
   characteristic in Borgomeo et al. (2015), Zaniolo et al. (2024), and
   the broader SSI literature (McKee et al. 1993; Vicente-Serrano et al.
   2010).
2. **Mean event cumulative deficit** (``mean_magnitude``) — mean over
   events of the cumulative SSI excursion (volume axis). Continuous;
   standard in severity-duration-frequency analysis (Yevjevich 1967).
3. **Time-in-drought fraction** (``time_in_drought_fraction``) — months
   with SSI ≤ −1 divided by total months (persistence axis). Continuous
   in [0, 1] and not subject to the integer-month clustering that
   disqualifies drought duration as a coverage axis.

Drought duration and peak-severity month are excluded because both
cluster at discrete monthly values and produce poor space-filling
properties on the resulting Pareto archives. The legacy three-tuple
(``mean_duration``, ``mean_avg_severity``, ``peak_severity_month``)
remains available as the ``legacy`` preset in
:data:`src.drought_metrics.PRESETS` for reproducibility of pre-2026-04-27
runs.

The accumulation period and the anti-ideal placement protocol use the
defaults inherited from the existing pipeline: SSI-3 monthly accumulation
and ``HEADROOM_TIMES_MAX`` placement at ``1.5 × historical_max`` for
unbounded metrics, ``CONSTANT = 1.0`` for the time-in-drought fraction.
DD-11 requires ``D_j ≤ D*_j`` for every feasible ``D``; the
:class:`AntiIdealRule.CONSTANT` rule on the time-in-drought metric
guarantees this trivially because the fraction is bounded above by 1.

**Alternative metric sets** are shipped in the registry for ablation:
``extreme_event`` (worst-event variants of severity and magnitude),
``trace_fdc`` (Q10 flow + time-in-drought + magnitude). Swapping is a
single ``--metric-set <name>`` CLI flag; see
:data:`src.drought_metrics.PRESETS`.

**Pointers:** :mod:`src.drought_metrics` (``DroughtMetric`` dataclass,
``REGISTRY``, ``PRESETS``, ``resolve_metric_set``, ``compute_anti_ideal``);
:mod:`src.experiment_config` (``metric_set`` field, default ``"primary"``);
:func:`src.objectives.compute_ssi_drought_characteristics` (extended with
``time_in_drought_fraction`` and ``q10_flow_neg`` keys).

---

## DD-05: Plausibility constraint family

*Status: SUPERSEDED by DD-14.*

The 2026-04-12 entry weighed several flow-statistic plausibility constraints
(lag-1 autocorrelation, non-drought mean, seasonal cycle, Hurst, cross-site
correlation) and settled on a soft-constraint regime under Deb (2000)
constraint domination. That regime is no longer the production choice.
Production constraint selection is in DD-14. The hydrologic five-statistic
regime is retained as one ablation arm in the SI.

---

## DD-06: SynHydro integration

*Status: SETTLED 2026-04-13.*

Production uses a thin wrapper (`src/kirsch_wrapper.py`,
`KirschBorgWrapper`) around SynHydro's validated Kirsch generator rather
than a standalone reimplementation. The wrapper exposes two DV injection
modes — an index mode that maps DVs to historical year indices, and a
residual mode that maps DVs to empirical-CDF residuals — and feeds either
through SynHydro's normal-score transform, Cholesky correlation step,
inverse transform, and Dec–Jan handling.

The choice between the two injection modes for production is an open
sub-decision tracked in `code_alignment_backlog.md` and exercised by
the wrapper-mode ablation experiments in `workflows/0N_<stage>/15_*` and
`16_*`.

**Pointer:** `src/kirsch_wrapper.py`.

---

## DD-07: Borg interface

*Status: SETTLED 2026-04-15. Updated 2026-04-28 -- MM Borg only.*

MM Borg MOEA (multi-master, MPI) is the **only** optimizer used anywhere
in the project. The analytic benchmark, the dimension sweep, the
epsilon x NFE sensitivity grid, the calibration runs, the ablations, and
the production Cannonsville runs all dispatch through
`src.borg_runner.run_optimization` and run via mpirun. The serial Borg
backend is retained as a single-process fallback for one-off interactive
debugging only. The earlier EpsNSGAII (platypus) stand-in is removed
from the codebase entirely as of the 2026-04-28 cleanup -- it was only
ever a laptop convenience for runs that should always have been done on
HPC. NSGA-III and MOEA/D remain out of scope by user direction.

**Operational topology rule.** The MM Borg wrapper picks
``n_islands`` automatically from the MPI rank count. Single-island
master-slave mode (``n_islands == 1``) intermittently SIGSEGVs on a
worker rank during MPI startup once the per-master worker count
exceeds ~6, so :func:`src.borg_runner._auto_islands` enforces
``n_islands >= 2`` for any allocation of >= 5 ranks. The production
stage 04 slurm matches the Hadka & Reed 2015 ratio (4 islands at 120
ranks). Per-call ``maxEvaluations`` is divided by ``n_islands`` inside
the wrapper because Borg's ``solveMPI`` argument is per-island, not
total. See ``docs/borg_integration_notes.md`` for the complete bug
analysis.

**Pointer:** `src/borg_runner.py`; DD-11 verification.

---

## DD-08: Zero-drought solutions on the Pareto front

*Status: DEFERRED.*

In the residual injection mode the Pareto front contains a single solution
at $D = 0$ corresponding to a trace with no detected drought events. This
is the mathematically expected corner of the feasible image and is not
pruned. If exclusion becomes desirable downstream, the recommended fix is
a hard constraint requiring at least one detected drought event per trace.
Not a priority at current stage.

---

## DD-09: Coverage baseline (library subsampling in drought characteristic space)

*Status: SETTLED 2026-04-13.*

The apples-to-apples baseline for the coverage comparison is a 10 k+-trace
Kirsch-Nowak library subsampled in drought characteristic space via LHS or
Sobol, not a QMC sample on the bounding box of drought characteristics. Only
the library baseline shares MOEA-FIND's feasible region; the bounding-box
QMC alternative samples points the generator cannot reach.

The closest published precedent is Bonham et al. (2024), who subsample in
generator *input* space rather than drought *characteristic* space; that
contrast is one of the paper's positioning points.

**Pointer:** `src/library.py`; `workflows/0N_<stage>/05_kirsch_library_build.py`;
`workflows/0N_<stage>/06_library_subsample_baseline.py`.

---

## DD-10: Feasible-region framing

*Status: SETTLED.*

MOEA-FIND's contribution is structured coverage of the *feasible* drought
hazard region, not raw uniformity over an unconstrained box. Coverage metrics
that compare MOEA-FIND to a QMC sampler on the bounding box are misleading
because the QMC sampler places mass in regions the generator cannot reach.
DD-09 fixes the comparison by restricting baselines to the same feasible
region. The paper must lead with the feasible-region framing rather than
with raw discrepancy values, or readers will misread the coverage tables.

Empirical evidence that the construction tiles the feasible region rather
than its boundary is in DD-11.

---

## DD-11: L1 formulation and verified interior coverage

*Status: SETTLED 2026-04-14 on the analytic benchmark; the high-dimensional
extension is the open DD-12.*

**Formulation.**

$$
f_j(x) = D_j(x) \quad \text{for } j = 1, \ldots, K, \qquad
f_{K+1}(x) = \lVert D(x) - D^* \rVert_1.
$$

Under the assumption $D_j(x) \leq D^*_j$ for every feasible $x$ and every
$j$, the absolute value reduces to $D^*_j - D_j$, the sum
$\sum_{j=1}^{K+1} f_j$ is the constant $\sum_j D^*_j$, every feasible
objective vector lies on a codimension-one affine subspace $S$ of the
$K+1$-objective space, and every pair of feasible points is mutually
non-dominated. Borg's epsilon-dominance archive tiles $S$, and the
projection onto the first $K$ coordinates is a bijection onto the feasible
drought hazard region in drought characteristic space.

The implementation in `src/objectives.py` uses this form (not the degenerate
$f_j = \lvert D_j - D^*_j \rvert$ form, which would collapse the front to
the anti-ideal). The affine identity is verified to machine precision on
the analytic benchmark; per-archive-member residuals are bounded by the
IEEE 754 rounding budget.

**Empirical interior-coverage verification.** Documented in
`evidence/shell_vs_interior_diagnostic.md`. A dimension sweep on a constrained
K-dimensional analytic benchmark shows that the archive matches reference
QMC samplers on mean L1 distance from the anti-ideal and on signed orthant
occupancy across the dimensionalities tested, and matches or exceeds them
on interior mass fraction. All numerical values are produced by MM Borg
MOEA via MPI on HPC (per DD-07).

**Anti-ideal placement.** The construction depends on $D_j(x) \leq D^*_j$
for every feasible $x$; if any archive member violates this in any
dimension, the affine identity fails silently. For the analytic benchmark
$D^*$ is placed strictly outside the feasible region by construction. For
the hydrologic case the anti-ideal placement protocol is the open
sub-decision flagged in `code_alignment_backlog.md` Item 8.

**Open extensions:**

- The high-dimensional Cannonsville case is DD-12.
- Behaviour on a non-convex feasible region better matching the real Kirsch
  case has not been tested in this DD; not blocking.

---

## DD-12: High-dimensional empirical audit on Cannonsville

*Status: PENDING. Hard gate on §3.2 prose.*

DD-11 verifies the construction on a low-dimensional convex analytic
benchmark. The Cannonsville case has hundreds of continuous decision
variables and a non-convex feasible drought hazard region produced by the
SSI event extraction pipeline (or by whatever drought-characterisation
pipeline DD-04 settles on). DD-11 does not substitute for verifying interior
coverage on the actual archive; the two diagnostics stress different limbs
of the construction.

**Failure modes the audit must rule out.**

1. *Shell-only coverage.* Archive concentrates at the outer boundary of the
   feasible region.
2. *Orthant collapse.* Archive misses one or more of the $2^K$ signed
   orthants around the anti-ideal.
3. *Sub-grid clustering.* Fine-scale clustering interior to the feasible
   region (the only mild deviation observed at the highest $K$ tested in
   DD-11; should not appear in the production Cannonsville case at the
   planned NFE budget).
4. *Disconnected feasible components.* Components reachable only through
   low-probability regions of the decision space, which the adaptive
   variational operators may or may not discover.

**Audit protocol.** Run after Phase β completes against a feasible-region-
restricted Latin hypercube subsample of the 10 k+-trace Kirsch-Nowak library
on the same drought characteristic axes. Compare distributions of L1 distance
from the anti-ideal, interior mass fraction, signed orthant occupancy,
nearest-neighbour CV, and feasible-region-restricted L2 star discrepancy.
Visual diagnostic by pairwise projection onto the planes of the drought
hazard space.

**Decision rule.** Audit passes if the archive matches or exceeds the
library-restricted reference on interior mass fraction and orthant occupancy
and if no disconnected feasible component is visible. Audit failure forces
§3.2 of the manuscript to report the failure mode and §4 to add a specific
caveat naming the dimensionality and epsilon budget at which it appears.

**Implementation.** Post-processing script analogous to
`workflows/0N_<stage>/diag_shell_vs_interior.py`, applied to the Phase β
archive and the Phase γ library.

---

## DD-13: Manuscript outline and prose rules

*Status: SUPERSEDED 2026-04-27.*

The 2026-04-14 outline (five top-level sections with a fixed seven-figure
sequence) is no longer authoritative. The current manuscript scope and
section structure live in `drafts/manuscript_main_draft.md` (now a
scaffold). The hard prose rules from this DD are in
`governance/style_guide.md` §§5.11–5.17 and `§6` and govern every edit to
the draft. Use those, not this entry.

---

## DD-14: Anderson-Darling DV-uniformity constraint

*Status: SETTLED 2026-04-21.*

Production runs use a single Anderson-Darling goodness-of-fit statistic on
the flattened decision-variable vector against $U[0, 1]^d$, with the
tolerance bootstrap-calibrated against uniform samples
(`workflows/0N_<stage>/diag_dv_uniformity_calibration.py`). The constraint
is a soft constraint under Deb (2000) constraint domination.

**Rationale.**

- *Calibration traceability.* A single scalar threshold with a clear
  operational definition: the constraint passes any DV configuration
  reachable from a uniform distribution at the calibrated significance
  level. Replaces the per-flow-statistic tolerance calibration the
  hydrologic regime required.
- *Tighter feasible region.* The L2-star alternative tolerated DV
  configurations whose generated traces saturated the cyclic peak-month
  metric at its upper bound, which is a numerical artefact rather than a
  physically meaningful drought. The AD statistic, more sensitive to
  distributional tails, eliminates this.
- *Comparable coverage.* The drought-characteristic coverage achieved
  under AD was within seed-level noise of the hydrologic baseline at
  production NFE, on the empirical comparison documented in
  `workflows/04_moea_find_single_site/dv_uniformity_ablation.py` and
  `14_dv_uniformity_compare.py`.

The hydrologic five-statistic regime and the L2-star DV regime are retained
as ablation arms in the SI.

**Caveat.** The empirical comparison underlying the rationale was conducted
with a particular set of drought metrics. The decision itself (use AD on
DV space) is independent of the drought metric choice in DD-04, but the
specific coverage numbers tied to that comparison are not preserved here
because DD-04 may change.

**Pointer:** `src/constraints_dv.py`;
`workflows/04_moea_find_single_site/run_moea_find.py` (production driver,
default `--constraint-mode dv_uniform --statistic ad`);
`workflows/0N_<stage>/diag_dv_uniformity_calibration.py`.

---

## DD-15: Joint metric-and-T justification protocol

*Status: PENDING (2026-04-29) — Stage 1 complete, Stage 2/3 in flight.*

DD-01 (trace length) and DD-04 (drought metrics) interact: short T
degrades event-derived metrics (e.g. SSI-12 events become rare in
5-yr blocks); short T may also alter the inter-metric correlation
skeleton. The current `primary` preset (DD-04) and `n_years_out = 20`
default (DD-01) were each settled by analogy to prior art; neither
choice rests on a numerical comparison across the joint (K, T) space.
The DD-15 protocol resolves both decisions with a single empirical
sweep before the next production MM Borg commit.

**Coarse-bracket T grid:** `T ∈ {5, 10, 20, 30}` water years.
Brackets the existing T=20 default with a hazard-period region (5, 10)
that the user prefers and a long-trace anchor (30).

**K cardinality scope:** K ∈ {3, 4} evaluated jointly with T. K=3
matches the production `primary` preset; K=4 is included to test
whether an additional axis adds enough hazard coverage to justify the
extra Borg dimension.

**Workflow location:** overhauled `workflows/02_calibration/` in place
(no directory split). Five chained stages, each its own SLURM script:

1. **Stage 1 — historical T-blocks.** For each T, run the 28-metric
   candidate library (`src/extended_drought_metrics.py`) on every
   stride-1 historical block, with SSI-3, SSI-12 and Q80 calibrators
   fit *once* on the full 73-WY record (DD-11 lock-in carries through).
   Driver: `workflows/02_calibration/t_sensitivity_historical.py`;
   aggregator: `t_sensitivity_aggregate.py`. Outputs include
   `degeneracy.csv` (the small-T viability headline),
   `cluster_stability.csv` (Adjusted Rand Index between cluster
   memberships across T-pairs), and `kset_stability.csv` (Jaccard of
   top K-sets vs T=20 reference).
2. **Stage 2 — Kirsch-Nowak fidelity.** Generate 10,000 baseline Kirsch
   realizations at every T (driver:
   `workflows/03_kirsch_library/build_library_extended.py`) and
   compare metric distributions to historical T-blocks via per-metric
   Kolmogorov–Smirnov statistics, Frobenius shift of the Spearman
   correlation matrix, and L2-star/NN-CV coverage in each candidate
   K-set's subspace
   (`workflows/02_calibration/t_sensitivity_kirsch_compare.py`).
3. **Stage 3 — joint K × T decision matrix.** For each of the top-5
   K=3 and K=4 strict-rung K-sets at each surviving T, score on a
   five-component composite (min robust spread, constraint-rung
   penalty, 1 − max |ρ|, 1 − KS_max, 1 − norm cost). Cost component is
   driven by per-T MOEA inner-loop timing
   (`workflows/02_calibration/eval_cost_timing.py`); decision matrix
   driver `workflows/02_calibration/decision_matrix.py` writes
   `pareto_front_KxT.json` with the recommended (K*, T*) and two
   alternates.
4. **Stage 4 — confirmatory MM Borg run at (K*, T*).** Production
   driver `workflows/04_moea_find_single_site/run_moea_find.py` with
   `--n-years T*` and a new `--metric-set` preset registered in
   `src.drought_metrics.PRESETS`. 120-rank, 4-island, NFE = 200,000.
   User reviews and approves the (K*, T*) recommendation before this
   stage commits compute.
5. **Stage 5 — manuscript figures.** Four figures
   (`src/plotting/02_calibration/figure_{a,b,c,d}_*.py`):
   metric stability across T (violins + Kirsch ridges + MOEA scatter),
   Spearman clustermap preservation across {Hist, Kirsch, MOEA},
   drought-space coverage with L2-star / NN-CV bar, and exemplar
   MOEA-FIND traces with SSI shading and Hist-IQR annotation.

**DD-11 carry-through.** All stages preserve the L1-Device guardrail:
`f_j = D_j + ‖D − D*‖_1`, never `|D_j − D*_j|`. Anti-ideal placement
remains `1.5 × historical max` for non-cyclic axes (none of the K=3/K=4
candidates considered are cyclic; legacy cyclic peak-month is excluded
by saturation diagnostics, consistent with DD-04). The new K* preset
inherits the same `AntiIdealRule.HEADROOM_TIMES_MAX` placement.

**Verification gates.** Three explicit decision points:

1. After Stage 1: Drop any T whose strict clusters-and-concepts K=3
   rung is empty OR whose hard-degeneracy fraction (zero/NaN > 25%
   per metric) exceeds 50% of all 28 candidates. Saturation alone
   does not drop a T (the spread-screen IQR>0 check already filters
   degenerate metrics individually).
2. After Stage 2: Drop any T whose KS_max > 0.10 or whose
   Frobenius correlation shift > 1.5.
3. After Stage 3: User reviews `pareto_front_KxT.json` before Stage 4.

**Carve-outs.** Multi-site Cannonsville-plus-Pepacton extension
remains deferred (DD-02 / publication plan §4); DD-15 covers the
single-site DRB Cannonsville production case study only. The
event-level framing (single-event traces, 1–8 yr) noted in DD-01
remains a follow-up direction and is *not* part of DD-15's T grid.

**Pointer:** `src/metric_screening.py`;
`workflows/02_calibration/{t_sensitivity_historical,t_sensitivity_aggregate,t_sensitivity_kirsch_compare,eval_cost_timing,decision_matrix}.py`;
`workflows/03_kirsch_library/build_library_extended.py`;
`src/plotting/02_calibration/figure_{a,b,c,d}_*.py`;
`workflows/02_calibration/run_t_sensitivity.sh` (orchestrator
blueprint); `outputs/02_calibration/decision_matrix/pareto_front_KxT.json`
(authoritative recommendation, written by Stage 3).

---

## DD-16: Magnitude-varying sensitivity in drought-hazard space

**Status.** COMPLETE (Stages A + B). Method, code, and tests landed
2026-04-30. K=3 archive proof-of-concept run completed 2026-04-30
(SLURM job 217983, 8 min wall, N=3308 full archive). Figures
generated 2026-04-30 (job 217984). Re-run on the post-DD-15 K*
archive is a one-line invocation once DD-15 Stage 4 completes.

**Decision.** Section 3.4 of the manuscript (and SI-12b) reports a
magnitude-varying sensitivity analysis (MV-SA) adapted from
Hadjimichael et al. (2020). The method is fixed as follows:

1. **Factor space** = the optimised MOEA-FIND objective axes
   (drought-hazard characteristics from the upstream archive's
   `objective_keys`, identical to Stage 08). Non-optimised
   characteristics are excluded because they have no space-filling
   guarantee. A uniform-random *control* factor is appended so the
   magnitude-varying noise floor is empirically visible.

2. **Magnitude axis** = a single Pywr-DRB outcome from the
   Stage-06 metric bank. Headline axis: `nyc_min_storage_frac`
   (NYC minimum combined reservoir storage fraction). Supporting SI
   axes: `nyc_drawdown_days_below_0.25` and
   `montague_flow_vulnerability`. Multiple axes are run as
   independent invocations; results never fuse outcomes onto a
   single sweep.

3. **Percentile grid** = 19 evenly-spaced points
   `{0.05, 0.10, …, 0.95}`. The 0% and 100% endpoints are
   excluded because the indicator collapses to a constant there.

4. **Response form** = within-trace percentile (Hadjimichael Variant
   2). For each trace *i* and percentile τ, `Y_i(τ)` is the τ-th
   percentile of trace *i*'s own 20-year annual NYC min storage
   series — continuous, no binary indicator, no external historical
   reference threshold. The 20 calendar years of each trace provide
   the inner ensemble that Hadjimichael obtains from stochastic
   realizations. Requires a precomputed
   `annual_nyc_min_storage_frac.parquet` (20 years × N traces) from
   Stage 09 `precompute_trace_series`. A *conditional* response on a
   secondary outcome (FFMP Level 6; Montague vulnerability) within a
   window of `window_frac=0.30` of the rank-τ slice is planned for
   SI-12b. The exceedance indicator `I(τ) = 1{M ≤ M_quantile(τ)}`
   was evaluated but rejected: at N=3308 the SALib Delta noise floor
   (~0.13) swamps real signal on Bernoulli Y, and the response
   collapses because the same variable drives both the axis rank and
   the response.

5. **Methods** = Delta moment-independent (manuscript anchor and only
   method run in the K=3 proof-of-concept). PAWN density-based and
   RBD-FAST are implemented in the engine but deferred to the K*
   re-run to contain runtime. Cross-method triangulation will be
   reported in SI-12b.

6. **Bootstrap** = 50 replicates for the proof-of-concept run; will
   increase to 200 for the K* production run. Seed 42; percentile
   CIs on headline index and rank-Spearman stability per slice.

**Why the factor / axis flip.** Hadjimichael et al. (2020) apply
MV-SA over uncertain *parameters* (climate perturbations, demand
multipliers, infrastructure variables) with a stakeholder-relevant
output magnitude. We exchange the factor space for the
drought-hazard characteristic space and the magnitude axis for a
Pywr-DRB operational outcome. The diagnostic question becomes
"which features of drought drive operational stress at each
severity?" rather than "which knobs drive shortage at each
severity?". Drought characteristics are observable to planners,
they are exactly the dimensions MOEA-FIND structures coverage over,
and a flipped MV-SA is therefore the natural validator of the
coverage objective: it shows whether the dimensions we elevated as
optimisation targets are the dimensions that empirically drive
downstream operational stress, and whether their relative
importance shifts with severity.

**Why MV-SA in addition to §3.3 scenario discovery.** The §3.3
classifier (gradient-boosted tree on Pywr-DRB FFMP Level-6 / Hashimoto
labels) answers a separability question — is operational failure
separable in the D-characteristic space? MV-SA answers a complementary
dimension-importance question — *which* D dimensions drive
operational stress *at which severity*? The two diagnostics share
the factor space and the archive but ask different questions; both
appear in §3.

**Reviewer-defence implications.** MV-SA in characteristic space is
a quantitative response to HC-3 (Critique 6, "drought is emergent")
and sharpens the Critique-14 novelty claim — a direct comparison
against Hadjimichael 2020 in which the methodological extension is
precisely the reframing of the factor space.

**K=3 proof-of-concept results (2026-04-30).** N=3308 full archive,
Delta only, 19 percentiles, n_bootstrap=50, nr=25. Clear rank
crossover at τ ≈ 0.40: `mean_magnitude` leads at τ ≤ 0.35 (severe
storage stress), `time_in_drought_fraction` leads at τ ≥ 0.40
(moderate/favorable). `mean_severity` consistent rank 2–3. Control
factor sits at Δ ≈ 0.05–0.09 throughout (well below real factors).
Wall time: 8 min (19 MPI ranks, single-node exclusive).

**Pointer:** `src/magnitude_varying_sa.py`;
`src/plotting/magnitude_varying_sa.py`;
`workflows/09_magnitude_varying_sa/run_mv_sa.py`;
`workflows/09_magnitude_varying_sa/configs/delta_within_trace.yaml`;
`workflows/09_magnitude_varying_sa/slurm/run_mv_sa.slurm`;
`tests/test_magnitude_varying_sa.py`.

---

## DD-15a: Tier-G hazard-clean candidate metrics (DD-15 follow-up)

*Status: PENDING (2026-04-29) — diagnostic re-run in flight at T = 5.*

User review of the first DD-15 (K*, T*) recommendation surfaced four
methodological complaints against the legacy 28-metric candidate pool.
The Tier-G additions in `src/extended_drought_metrics.py` resolve the
three valid complaints; the fourth is verified as already-correct.

**Complaint 1 — SSI calibrator must be fitted on the full record.**
Verified already correct: `t_sensitivity_historical.py:160-165` and
`build_library_extended.py:80-92` fit SSI-3 and SSI-12 once on the
full 73-WY Cannonsville record and `.transform()` per T-block; Q80 is
also fixed at the full-record 20th-percentile. No code change needed
(carries DD-11 lock-in through Tier G unchanged).

**Complaint 2 — `cv_annual_min` is block-normalised and uninformative
for cross-block comparison.** Confirmed: the metric divides
`std(annual_min_in_block)` by `mean(annual_min_in_block)`, so a wet
and a dry block can collide at the same CV. **Fix:** drop
`cv_annual_min` from the K-set enumeration pool; replace with
`min_annual_zscore`:

```
min_annual_zscore = (refs.annual_mean_mean - block_min_annual_mean)
                  / refs.annual_mean_std
```

where `refs` is :class:`src.extended_drought_metrics.FullRecordRefs`,
constructed once per pipeline invocation via
`FullRecordRefs.from_full_record(monthly_1d_full)`. Larger value =
block contains a much-drier-than-typical year.

**Complaint 3 — mean-event metrics obfuscate extremes.** A 5-yr block
with one extreme plus one minor event lands `mean_severity` halfway
between, losing operational interpretation. **Fix:** add two integral
/ rolling-window replacements that don't extract events:

* `total_deficit_ssi3 = -Σ(SSI-3 over deficit-months)` — every
  deficit month contributes its full depth additively. No event
  extraction. Continuous, full-record-anchored via SSI calibration.
* `min_36mo_ssi3 = -min(rolling_36mo_mean(SSI-3))` over the block —
  captures the single worst sustained sub-period without picking just
  the worst single month and without saturating the way
  `worst_severity` does at small T.

**Complaint 4 — `mean_recovery_time` is integer-stratified.**
Confirmed: recovery time is in months (integer). Same DD-04
exclusion logic that dropped `mean_duration` and
`peak_severity_month` for clustering at integer values. **Fix:** drop
from the enumeration pool; the new sustained-severity and integrated-
deficit metrics cover the recovery / persistence axis continuously.

**The Tier-G `HAZARD_CLEAN_METRICS` allowlist** for K-set enumeration:

| Metric | Concept | Definition |
|---|---|---|
| `total_deficit_ssi3` | volume_integral | −Σ(SSI-3) over deficit months |
| `min_36mo_ssi3` | sustained_severity | block-min of 36-month rolling SSI-3 mean (negated) |
| `min_annual_zscore` | annual_min | z-score of driest annual mean against full record |
| `q10_zscore` | flow_tail | z-score of block q10 against full-record monthly |
| `q25_zscore` | flow_tail | z-score of block q25 against full-record monthly |
| `q80_total_deficit` | volume_integral_q80 | Σ(Q80_full − flow) over flow < Q80_full months |
| `time_in_drought_fraction` | persistence_share | months SSI-3 ≤ −1 / total months (kept) |

All seven metrics are continuous, full-record-anchored, non-averaged
across events, and not stratified at integer values. Sign convention:
larger value = more drought-stressed (matches DD-11 L1 device).

**Decision-matrix selection mode.** `decision_matrix.py` gains a
`--metric-pool {full,hazard-clean}` flag (default `hazard-clean`).
The hazard-clean pool restricts K-set enumeration to the seven
metrics above; the full pool retains the existing 34 candidates
subject to the saturation/zero-degeneracy/concept filters introduced
in DD-15. The legacy pool remains visible in the DF1 / DF2
diagnostic figures so the reader can see *why* each excluded metric
was excluded.

**No regeneration required to re-evaluate.**
`workflows/03_kirsch_library/recharacterize_library.py` reads the
existing `library.npy` per T and re-runs `compute_all_candidates`
with the updated metric set, writing a refreshed
`characteristics_extended.npz`. ~25 min per T at N = 10 000; avoids
the ~30–90 min Kirsch generation step.

**T = 5 user-locked.** The diagnostic re-run only re-evaluates K-set
candidates at T = 5; the T = 10 / 20 / 30 historical and Kirsch
artefacts remain on disk and feed the SI sensitivity-to-T section
unchanged.

**Pointer:** `src/extended_drought_metrics.py:FullRecordRefs`,
`compute_hazard_clean_metrics`, `HAZARD_CLEAN_METRICS`;
`workflows/03_kirsch_library/recharacterize_library.py`;
`workflows/02_calibration/decision_matrix.py:--metric-pool`;
`figures/02_calibration/diagnostics_for_decision/df{1,2,3,4}_*.png`.

---

## DD-15b: Short-block (T = 1, T = 2) reframe under investigation

*Status: PENDING (2026-04-29) — diagnostic suite running; manuscript
reframe deferred until the data supports the move.*

User review of the DD-15a hazard-clean K=3 / K=4 recommendation
exposed two limits of the multi-year framing for the DRB Cannonsville
case study:

1. **Even the hazard-clean metrics rank-correlate at median |ρ| = 0.65
   (Kirsch) / 0.82 (Hist) at T = 5.** The historical 5-yr blocks
   project onto a near-1-D axis — variation across blocks is dominated
   by "is this block in or near the 1960s drought?". At cap = 0.6,
   strict-rung K = 3 admits zero candidates; K = 4 is far worse.
2. **The DRB physical system is not multi-year-persistence-limited.**
   NYC's combined operating capacity is small relative to annual
   demand; sub-annual stress patterns (low summer recession, dry July–
   August) drive operator vulnerability. The literature default of
   T = decades does not match the physical scale of the system.

The **DD-15b reframe** is a candidate course correction: shift the
sampling unit from "structured multi-year drought" to "structured
stress year (or stress 2-year period)", with antecedent factors
(initial storage, demand, prior-year hydrologic class) treated as
**orthogonal post-hoc LHS factors** that combine factorially with
the MOEA-FIND hydrologic samples in the §3.3 / §3.4 SA arm.

**T grid** for the diagnostic: `T ∈ {1, 2}`, with the existing T = 5
hazard-clean results retained as a reviewer-comparison anchor.

**Metric pool** (`src/short_block_metrics.py:SHORT_BLOCK_METRIC_NAMES`):
12 candidates spanning two tiers:

* Tier H (raw flow, no SSI fitting): `total_flow_neg`, `jja_total_neg`,
  `djf_total_neg`, `min_monthly_flow_neg`, `min_3mo_rolling_neg`,
  `summer_recession`, `min_annual_zscore`, `q10_zscore`, `q25_zscore`.
* Tier I (SSI-3 single-block with 3-month burn-in):
  `min_ssi3_neg`, `total_deficit_ssi3_block`,
  `time_in_drought_fraction_block`.

All 12 are continuous, full-record-anchored where reference-fitting is
required (SSI-3 calibration uses the full 73 WY record), non-averaged
across events, and not stratified at integer values. Larger value =
more drought-stressed (consistent with DD-11 L1 device).

**SSI-3 burn-in protocol.** SSI-3 monthly values require three months
of preceding flow data. The diagnostic prepends three burn-in months
to every evaluation block (historical or synthetic) and computes
SSI-3 across the burn-in + evaluation window; metric aggregation uses
only the evaluation window. For the production MOEA driver, this
implies a 15-month DV vector at T = 1 (3 burn-in + 12 evaluation) and
a 27-month DV vector at T = 2.

**Antecedent architecture.** Antecedents (initial storage, demand,
prior-year hydrologic class) are **not** Borg DVs and **not** part of
the hazard space coverage objective. They enter the framework as a
post-hoc LHS factorial sample combined with the MOEA-FIND archive
during the Pywr-DRB simulation arm (§3.3 / §3.4). This separates the
"what hydrology stresses this system" question from the "given a
stressing hydrology, what initial conditions and demands matter"
question — the former is Borg-search, the latter is variance-based
sensitivity analysis.

**Headline result so far** (from historical-only correlations,
pending Kirsch confirmation):

| Setting | n metrics | median \|ρ\| | strict-rung K=3 @ cap = 0.6 | K=4 @ cap = 0.6 |
|---|---:|---:|---:|---:|
| T = 1 short-block | 12 | 0.713 | **10** | 1 |
| T = 2 short-block | 12 | **0.550** | **66** | **56** |
| T = 5 hazard-clean | 7 | 0.817 | 0 | 0 |

T = 2 is provisionally the sweet spot: the metric pool admits dozens
of viable K = 4 strict-rung sets at cap = 0.6, vs zero at T = 5. T = 1
has fewer options because `total_flow_neg` and `min_annual_zscore` are
deterministic transforms of each other (perfect correlation).

**No manuscript reframe yet.** Per user request, the framing shift is
deferred until the Kirsch ensemble confirms the historical pattern
and the joint K-set / cost / coverage diagnostic supports the move
end-to-end.

**Pointer:** `src/short_block_metrics.py`,
`workflows/02_calibration/short_block_screening.py`,
`workflows/03_kirsch_library/build_short_block_library.py`,
`workflows/02_calibration/short_block_kset_table.py`,
`src/plotting/02_calibration/diagnostics_short_block.py` (DF7–DF9);
`outputs/02_calibration/short_block_screening/T{1,2}/`,
`outputs/03_kirsch_library/build_short_block_library/n10000_s42/`.

---

## DD-15c: T=1 first-run reveals metric-formulation pathologies

*Status: PENDING (2026-04-30) — T=1 reframe is retained; the K=4 metric
set chosen for the first production run produced symptoms inconsistent
with the DD-11 L1 device. The metric set will be reformulated; the
investigation below documents the observations and constrains the
reformulation. Not a reversal of DD-15b — T=1 stays.*

The first end-to-end T=1 short-block production runs (`stage 04`,
2026-04-30; jobs 217986, 217988–217990) used the K=4 set
`{djf_total_neg, summer_recession, aug_zscore, ond_total_neg}` from
`PRESETS["short_block_drb"]` in `src/drought_metrics.py`. A four-way
ablation (residual+AD baseline, residual+AD anti-ideal-iter2, residual+KS,
index+AD) at 200k NFE / 120 ranks each completed in 73–103 s per run
(T=1 raw-flow eval ≈ 40 ms vs ~300 ms at T=20-with-SSI).

The archives are well-formed by the DD-11 hyperplane diagnostic
(mean of $\Sigma_j D_j(x) + \|D-D^*\|_1$ within 0.6% of $\Sigma D^*_j$
in all four variants), but exhibit four interrelated pathologies that
trace back to the metric definitions:

### 1. Sign-asymmetry / unboundedness on the non-drought side

`djf_total_neg`, `ond_total_neg`, and `summer_recession` are computed as
*negated* raw-flow aggregates so that "larger value = more drought-
stressed" (DD-11 sign convention). The consequence is that the
**non-drought side of these metrics is unbounded** — a wet DJF can
produce `djf_total_neg = −10⁴` cfs·month, with no physical floor short
of arbitrarily large flow. The drought side is bounded (zero flow → 0
for `_total_neg` metrics; finite for `summer_recession`). The metric
range is therefore $(-\infty,\, 0]$ rather than the bounded
$[0, D^*_j]$ implicitly assumed by HEADROOM_TIMES_MAX anti-ideal
placement.

`aug_zscore` is doubly unbounded — both tails are unbounded in
principle, with the drought side bounded by zero-flow August
($z \approx \mu_{aug}/\sigma_{aug} \approx 1$), and the wet side
unbounded above (the most extreme Pareto member reached $z = -36\sigma$
under KS, corresponding to August flow ≈ 30× its historical mean).

### 2. Magnitude / unit heterogeneity

The four-metric L1 device aggregates objectives with mismatched scales:

| Metric                | Units      | Historical range          | iter1 Pareto range        |
|-----------------------|------------|---------------------------|---------------------------|
| `djf_total_neg`       | cfs·month  | [−10000, −805]            | [−5938, −306]             |
| `summer_recession`    | cfs/month  | [−1000, 549]              | [−9224, +2263]            |
| `aug_zscore`          | σ (dim'less) | [≈ 0, 0.76]              | [−24.6, 0.74]             |
| `ond_total_neg`       | cfs·month  | [−9000, −196]             | [−8407, −119]             |

The Manhattan norm $\|D - D^*\|_1$ is dominated by the cfs·month
metrics (10³–10⁴) and is essentially blind to the dimensionless
`aug_zscore`. Per-axis epsilons in `src/drought_metrics.py` (500, 50,
0.1, 500) attempt to compensate but cannot fix the L1-norm aggregation
itself.

### 3. Anti-ideal placement breaks under (1) + (2)

`HEADROOM_TIMES_MAX` (1.5 × historical max) places $D^*$ outside the
historical envelope but **inside** the feasible envelope under T=1
short-block exploration. The MOEA pushed `summer_recession` to a
Pareto-max of 2263 (174% above $D^* = 824$) in iter1; the protocol
re-anchor (iter2: $D^* = 1.5 × 2263 = 3394$) was violated again by
the same Pareto-max behavior (4147; 22% over). At least 2–3 more
iterations would be needed to converge $D^*$, and convergence is not
guaranteed if the metric is unboundedly explorable on the
drought-stressed side.

`djf_total_neg` and `ond_total_neg` were patched to use
`AntiIdealRule.CONSTANT` with `anti_ideal_constant=0.0` (zero-flow
physical bound) during implementation; this kept those two coordinates
of $D^*$ stable but did not address the unbounded *negative* tail of
the metrics, which is what the FIND framework explores when the MOEA
is asked to find "diverse extremes".

### 4. Pareto front includes flood-shaped non-drought scenarios

By construction, the FIND L1 device finds K corner-extreme scenarios
plus interior trade-offs. With the K=4 metric set above, the corner
that maximizes one metric typically *minimizes* others into the
non-drought direction. The "max `aug_zscore`" trace in iter1 has
$z = -24.6$, August flow = 6189 cfs (10× the historical August mean):
this is a flood-summer scenario that is extreme in the *non-drought*
direction of `aug_zscore`. Pareto-archive monthly flows reach 33,020
cfs vs historical max ≈ 7700, with 0.85% of monthly values exceeding
10,000 cfs.

The DV-uniformity constraint did not catch these because the
constraint operates in DV space, not flow space; index mode
(structurally bounded to historical residuals) reduced but did not
eliminate the issue (max 34,771 cfs, 0.47% > 10,000).

### Diagnosis (user, 2026-04-30): metric definitions are the cause

The four pathologies above co-arise when the metrics are (i) raw flow
aggregates rather than bounded distributional statistics, and (ii)
mixed in scale and unit. **The reformulation must keep T=1 but
replace the metric set with bounded, scale-comparable quantities.**

Constraints on the reformulation, agreed 2026-04-30:

1. **T=1 stays.** The DD-15b reframe (sub-annual stress patterns; NYC
   is not multi-year-persistence-limited) is retained as the framing.
2. **Metrics must be bounded** on both the drought and non-drought
   sides, so that $D^*$ is a true extremum and the Pareto front does
   not extend into flood-shaped corners.
3. **Metrics should be scale-comparable** so the L1 device aggregates
   them on equal footing (per-axis epsilons can fine-tune but should
   not load-bear).
4. **Metrics should remain sub-annual** and DRB-physical (DJF, summer,
   August, OND windows are appropriate; the lesson is in *how* they
   are summarized, not *which* months).

Candidate reformulations (not yet evaluated; each produces bounded
distributional statistics):

* **Empirical-CDF positions** in $[0, 1]$ against the full-record
  per-window distribution. E.g. `djf_cdf` = empirical CDF position of
  this trace's DJF total flow within the historical DJF-total
  distribution. Drought side ≈ 0, normal ≈ 0.5, wet ≈ 1. $D^*$ is
  the corner of the unit cube; the Pareto front lies inside it.
* **Z-scores clipped to $[\,−c,\, +c\,]$** for some envelope (e.g.
  $c = 4$). Preserves interpretability; truncates the
  flood-extreme tail.
* **Standardized deficit volumes** scaled by historical IQR or
  monthly $\sigma$ rather than left in raw cfs·month.

The choice is open and should be made before the next stage-04
production submission. The four-archive ablation set serves as the
"before" baseline for any reformulated set.

**Why:** The DD-11 L1 device (and FIND geometry generally) presumes
metrics with bounded feasible support and comparable scale. Raw-flow
aggregates do not satisfy either presumption at T=1, where a single
month's value can drive an entire metric.

**How to apply:** Until DD-15c is resolved, the T=1 archives at
`outputs/04_moea_find_single_site/run_moea_find/{residual,index}_T1_*/`
are diagnostic-quality only. Do not feed them into stage 06 (Pywr-DRB
re-eval), stage 07 (scenario discovery), or stage 09 (MV-SA) without
filtering by per-objective bounds.

**Pointer:** `src/drought_metrics.py:PRESETS["short_block_drb"]`,
`src/short_block_metrics.py:SHORT_BLOCK_METRIC_NAMES`;
`outputs/04_moea_find_single_site/run_moea_find/`
{`residual_T1_nfe200000_s42_constrained_cmdv_uniform_stad`,
`...sfxiter2_stad`, `...stks`,
`index_T1_nfe200000_s42_constrained_cmdv_uniform_stad`};
diagnostic figures at the matching paths under
`figures/04_moea_find_single_site/run_moea_find/<slug>/t1_archive_diagnostics.pdf`.
