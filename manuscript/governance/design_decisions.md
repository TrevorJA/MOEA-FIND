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

*Status: PARTIALLY SETTLED (2026-04-12; revisit after DD-04).*

The production case study is single-site Cannonsville with moderate-length
traces (on the order of decades), matching the trace-level framing of
Borgomeo et al. (2015), Zaniolo et al. (2024), and Wheeler et al. (2025).
Aggregate drought metrics computed across the full trace.

An event-level framing in which short traces (1–8 years) carry a single
designed drought event is a possible follow-up direction. Event-level
metrics (peak intensity, cumulative severity, recovery rate) are natural
substitutes for trace-level aggregates and would shrink the decision-variable
count substantially. Not the production focus.

The exact trace length and the trace-vs-event question are tied to the
choice of drought metrics in DD-04 and may shift with that decision.

**Pointer:** `drafts/manuscript_main_draft.md` §2.4;
`workflows/04_moea_find_single_site/run_moea_find.py`.

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

*Status: SETTLED 2026-04-27 (production default; user review remains
open before HPC commits).*

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
