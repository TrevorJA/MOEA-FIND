# MOEA-FIND verification criteria (draft v0.2, DD-11-aligned)

**Status**: draft, 2026-04-17. Iterate with the user.

## Formulation lock-in (2026-04-17)

The DD-11 L1 Device is the locked formulation for MOEA-FIND. Both
entry points
([src/objectives.py::analytic_objectives](../src/objectives.py) and
[src/objectives.py::drought_objectives](../src/objectives.py))
implement

    f_j(x)     = D_j(x)                  for j = 1..K
    f_{K+1}(x) = ||D(x) - D*||_1

with `D*` placed outside the feasible region
(`1.5 × max_hist` for non-cyclic metrics, `12 × 1.5 = 18` for cyclic
calendar-month metrics) so that `D_j ≤ D*_j` holds for every feasible
`D`. Under that assumption the hyperplane identity
`sum_{j=1..K+1} f_j = sum_j D*_j` is constant across the feasible set
and every feasible point is Pareto non-dominated.

Two independent verifications are on record:

1. **Analytic interior coverage (K=2..6 cube).** Via
   `workflows/diagnostics/diag_shell_vs_interior.py --feasible-shape cube`
   as a SLURM job array, MOEA-FIND matches uniform-in-cube on mean L1
   distance from `D*`, hits 100% of `2^K` signed orthants at every
   tested `K`, and shows only a mild budget-driven shell lean at
   K=5..6. Summary figure at
   [figures/figSI_shell_interior_sweep.pdf](../figures/figSI_shell_interior_sweep.pdf);
   per-K figures at `figures/figSI_shell_interior_k{K}.pdf`; table at
   [outputs/diag_shell_vs_interior/sweep_table.md](../outputs/diag_shell_vs_interior/sweep_table.md).
2. **Hydrology hyperplane numerical check.** The Cannonsville 20-k
   NFE run (job 195334) reports `std/mean = 2 × 10⁻³` on the
   `sum(f_j)` identity, confirming the DD-11 algebra is numerically
   stable under the Kirsch-Nowak evaluate callback and MM Borg.

Any code change that touches `drought_objectives`,
`compute_ssi_anti_ideal`, or the cyclic-month D* placement must
preserve DD-11 semantics — regenerate the analytic SI figures as the
regression test. The guardrail is also captured in
[memory/feedback_dd11_l1_device.md](../../../.claude/projects/-home-fs02-pmr82-0001-tja73-Research-DRB-Pywr-DRB-MOEA-FIND/memory/feedback_dd11_l1_device.md).

## Purpose

Define the bar MOEA-FIND must clear before it is considered a working
drought-scenario generator for the Delaware River Basin analysis. Each
criterion must have a pass/fail test tied to a concrete code location so
that a future run can be judged mechanically.

## Method reference (load-bearing)

The authoritative specification of the MOEA-FIND objective formulation is
``manuscript/design_decisions.md §DD-11``. In summary, Borg minimises the
``K+1``-vector

    f_j(x)     = D_j(x)                    for j = 1..K
    f_{K+1}(x) = ||D(x) - D*||_1           (Manhattan L1 distance)

with ``D*`` placed outside the feasible region
(``1.5 × max_hist`` for non-cyclic metrics, ``12 × 1.5 = 18`` for
cyclic-month metrics). Under the resulting ``D_j ≤ D*_j`` assumption,
``sum_{j=1}^{K+1} f_j`` is constant on the feasible set, so every
feasible point is Pareto non-dominated and Borg's ε-box archive tiles
drought-characteristic space uniformly. Implementation at
[src/objectives.py::drought_objectives](src/objectives.py); anti-ideal
placement at
[src/experiment_utils.py::compute_ssi_anti_ideal](src/experiment_utils.py).

*Reading the Pareto.* Under this formulation the Pareto archive spans
from the "no drought" corner (``D_j → 0``, wet traces) through the
"severe drought" corner (``D_j → D*_j``, dry traces). An **ensemble
median** aggregated across the full archive is therefore not the right
comparator — it averages the no-drought and severe-drought ends together
and lands near the middle, which for Cannonsville reads as "slightly
wetter than historical" because the low-drought end of the front is
dense with wet traces. Per-criterion tests below use **subset filters**
on drought metrics before comparing statistics to the historical block
envelope; see Criterion 2 for the canonical example.

## Baseline conventions

- **T** = synthetic trace length in water years. Currently configured
  default is 20, but this default is *under active reconsideration*
  via the DD-15 joint metric-and-T justification protocol
  (`manuscript/governance/design_decisions.md` §DD-15). The candidate
  T grid for that sweep is `{5, 10, 20, 30}`. Until DD-15 lands the
  recommended (K*, T*) tuple in
  `outputs/02_calibration/decision_matrix/pareto_front_KxT.json`, this
  doc continues to assume T = 20 for verification numbers.
- **N_hist** = length of historical record in water years. At Cannonsville
  USGS 01436000, [src/experiment_utils.py:59](src/experiment_utils.py#L59)
  reports 73 water years.
- **Historical block ensemble**: the set of all overlapping `T`-year
  windows of the historical record (stride 1 by default). At the current
  Cannonsville record, `N_hist - T + 1 = 54` blocks for `T=20`. This
  ensemble is the primary comparator for synthetic traces; any
  "historical FDC / ACF / seasonal cycle" line should be rendered as the
  *envelope* of these blocks, not a single window, so that the comparison
  is not biased by a particular block's drought content. Reference
  implementation: `src/historical_blocks.py` (TBD in Phase C of
  `plans/gleaming-noodling-pearl.md`).
- **Drought classification**: SSI-3, calibrated **once on the historical
  record** and applied via `.transform()` to every synthetic trace. The
  calibration happens in
  [src/experiment_utils.py:260-262](src/experiment_utils.py#L260-L262)
  (`run_experiment` pre-fits `prefitted_ssi` on `monthly_1d` before the
  MOEA loop) and the same object is reused in the evaluation callback
  ([src/experiment_utils.py:289](src/experiment_utils.py#L289)), the
  Pareto trace regeneration
  ([src/experiment_utils.py:392](src/experiment_utils.py#L392)), and all
  constraints that need SSI
  ([src/experiment_utils.py:298](src/experiment_utils.py#L298),
  [src/experiment_utils.py:398](src/experiment_utils.py#L398)). The
  Kirsch-library evaluator in
  [workflows/experiments/05_kirsch_library_build.py:87-89](workflows/experiments/05_kirsch_library_build.py#L87-L89)
  and the constraint-calibration diagnostic in
  [workflows/diagnostics/diag_constraint_calibration.py:356-358](workflows/diagnostics/diag_constraint_calibration.py#L356-L358)
  follow the same pattern. There are no callers of `compute_ssi(flows)`
  without `reference_flows=` on synthetic traces. Drought events are
  contiguous months with SSI-3 < 0 per SynHydro's
  `get_drought_metrics`; "critical" events require SSI ≤ -1.

## Criteria

### 1. Drought-space coverage

**Statement.** The Pareto archive's `drought_metrics` cloud must (a)
overlap the historical block ensemble's drought cloud, (b) extend beyond
it toward the anti-ideal in all minimised dimensions, and (c) **not**
collapse onto a single corner.

**Pass test.** Compute historical per-block `(mean_duration,
mean_avg_severity, peak_severity_month)` → historical cloud `H`. The
Pareto archive's `drought_metrics` must have:

- bbox overlap with H > 50 % of H's volume (sanity: we don't wander off)
- at least one Pareto point with each objective above the historical 95th
  percentile (progress toward the anti-ideal)
- each objective's Pareto-range ≥ 2 × H-range (the archive is more
  diverse than history, as expected from an anti-ideal optimiser)

**Where to implement.** New
`workflows/diagnostics/verify_drought_coverage.py`; inputs are
`results.json` + historical block ensemble. Not yet written.

### 2. Low-flow directionality

**Statement.** For the subset of Pareto traces whose `mean_duration` is
at or above the historical block median, the ensemble FDC at exceedance
probabilities ≥ 70 % must lie *at or below* the historical block-FDC
envelope. A drought-seeking MOEA must not systematically inflate
low-flow months.

**Pass test.** Define `drought_subset` = Pareto traces with
`mean_duration ≥ median(H.mean_duration)`. Compute the FDC of each trace
and of each historical block. For every exceedance probability in
{70 %, 80 %, 90 %, 95 %, 99 %}:

- synthetic 50th percentile FDC ≤ historical block-FDC 50th percentile **and**
- synthetic 10th percentile FDC ≤ historical block-FDC 10th percentile

**Current status.** FAIL: at 90 % exceedance the whole synthetic
ensemble is 2-3× above the single-window historical line (see
`outputs/exp04_kirsch_single_site/residual_T20_nfe2000000_s42_constrained/figures/fig06b_fdc.pdf`).
Note the historical bar here is a single 73-year FDC, not a block
envelope — see Criterion 4 for why that still doesn't rescue the run.

**Where to implement.** Extend the new drought-coverage diagnostic.
Reuse `src.plotting.trace_diagnostics.plot_flow_duration_curve`'s FDC
helper but drop the fixed `[0.01, 0.99]` interpolation grid
([src/plotting/trace_diagnostics.py:267](src/plotting/trace_diagnostics.py#L267))
so tail percentiles are honest.

### 3. Hydrologic-statistical fidelity at non-drought timescales

**Statement.** For the subset of Pareto traces whose `mean_duration` is
**below** the historical block median (the "nominal" half of the
archive), the following must lie inside the historical block-ensemble
envelope:

- FDC in the 0-30 % exceedance band (high-flow regime).
- Monthly lag-1 autocorrelation.
- Monthly seasonal cycle (mean and std per month).

**Pass test.** For each of the three statistics, ≥ 90 % of
`nominal_subset` traces must fall inside the historical block envelope
at each sampled index. Statistics with a scalar per trace (e.g.
lag-1 AC) must have the `nominal_subset` distribution's
[5 %, 95 %] interval overlap the historical block
[5 %, 95 %] interval.

**Where to implement.** Reuse existing
`src.constraints._lag1_ac` and the per-month helpers in
`src.plotting.trace_diagnostics`.

### 4. No constraint-artifact compression of flow space

**Statement.** The Pareto archive's distribution of annual mean flows
must not be materially narrower than the historical block ensemble's
distribution.

**Pass test.** Let `A_syn` = per-Pareto-solution annual mean (in the
convention of
[src/constraints.py:229](src/constraints.py#L229) — sum of 12 monthly
cfs, averaged across years). Let `A_hist` = per-historical-block annual
mean. Require:

- **Match the historical spread** (user answer, open question 3):
  `0.9 ≤ std(A_syn) / std(A_hist) ≤ 1.1`. A narrower spread implies the
  constraint is compressing flow space; a wider spread implies the
  objective is pulling solutions outside the plausible range.
- **Corner behaviour consistent with DD-11.** The low-drought corner of
  the Pareto (``mean_duration ≤ historical median``) should have
  ``A_syn`` distribution skewed *wetter* than historical (these are
  wet, drought-free traces), and the high-drought corner
  (``mean_duration ≥ historical median``) skewed *drier*. If both
  corners cluster at a single ``A_syn`` value, the constraint has
  collapsed flow-space diversity and Criterion 4 fails.

**Current status (2 M-NFE archived).** Unresolved. Diagnostic
[workflows/diagnostics/diag_constraint_pull.py](workflows/diagnostics/diag_constraint_pull.py)
reported the 10 lowest-duration solutions cluster at +15-18 % above
`hist_annual_mean`; that observation is consistent with DD-11 (low-
drought traces should be wetter than historical). The `corr(mean_duration,
annual_mean) = -0.26` sign is correct for drought-seeking. The test as
now written (match spread, bifurcation across corners) has not been run
against the archived Pareto yet and is the pending empirical check.

**Where to implement.** Extend
[workflows/diagnostics/diag_constraint_pull.py](workflows/diagnostics/diag_constraint_pull.py)
or fold into the new drought-coverage diagnostic.

### 5. Drought seasonality

**Statement.** Drought onset months in the synthetic ensemble should
track the historical block-resampled drought-onset distribution.
Cannonsville historical droughts concentrate in summer/early autumn;
the synthetic archive must not produce systematically different onset
seasonality.

**Pass test.** For each Pareto trace and each historical block,
compute a histogram of drought onset months (12 bins). Chi-squared
test between the `drought_subset` aggregate histogram and the
historical aggregate histogram must have `p ≥ 0.05`.

**Where to implement.** Not yet — new diagnostic module TBD.

### 6. Multi-site spatial coherence (deferred)

**Statement.** Only applicable once Script 08 (multi-site MOEA) runs.
Cross-site correlation of drought-onset dates within each synthetic
trace must match the historical inter-site lag structure.

**Where to implement.** TBD.

## 2026-04-17 session notes

### Was the archived 2 M-NFE FDC inversion actually a methodology failure?

Provisional read: **probably not** — the FDC comparison was scored
against the wrong baseline.

- SSI-3 calibration was never per-trace (empirically confirmed; the
  lowest-duration Pareto trace scored through a fresh historical-fit SSI
  returns ``n_events=0`` identical to the stored ``pareto_chars``).
- The objective formulation was and now is DD-11 (``f_j = D_j``, ``f_{K+1}
  = ||D - D*||_1``). Under this formulation the Pareto archive is
  *supposed* to span the "no drought" corner (wet traces, ``D_j ≈ 0``)
  through the "severe drought" corner (dry traces, ``D_j ≈ D*_j``). An
  ensemble median drawn across the full archive therefore lands wetter
  than historical — that is an aggregation artifact, not a failure of
  the method.
- A brief code detour between 2026-04-17 afternoon and that evening
  replaced DD-11 with the degenerate ``f_j = |D_j - D*_j|`` variant that
  DD-11 explicitly rejects. A login-node NFE=100 smoke test confirmed
  DD-11's prediction: the Pareto collapsed to a tight cluster near
  ``D*`` (``mean_duration ∈ [45, 47]``, ``mean_avg_severity ∈ [10.0,
  10.7]``), exactly the ``D = D*`` attractor DD-11 warns about. The
  detour has been fully reverted and the tests now lock DD-11 in.

### What still needs empirical verification

The 2 M-NFE archive is in
``outputs/archive/exp04_kirsch_single_site/residual_T20_nfe2000000_s42_constrained``.
Before declaring it good or bad we need to score it against the
DD-11-aligned criteria above:

- **Criterion 2** (drought_subset FDC ≤ historical block envelope at
  high exceedance) — expected to pass if the DD-11 Pareto is tiling
  correctly.
- **Criterion 4** (flow-space spread not compressed by constraints) —
  the `diag_fdc_inversion` scatter already showed a +15 % wet bias at
  the low-drought corner. Under DD-11 that is not automatically a
  problem (the low-drought corner should be wet), but whether the
  constraint band actually binds feasibility needs an asymmetric check.

Phase-D unblock depends on these tests.

## Open questions for the user

1. For the multi-block historical envelope (Criterion 2-5), is stride-1
   overlap the right default, or a larger stride (e.g. T/2) to reduce
   correlation between blocks? **A: Stride 1 is fine, and will give better convergence.**
2. Does the "nominal vs drought subset" split in Criterion 3 need to be
   specific to each objective (`mean_duration ≥ median_hist`,
   `mean_avg_severity ≥ median_hist`, `peak_severity_month ≥ median_hist`
   — intersection or union)? **A: I'll decide from justifiability.**
   *Working rule*: use the **union** for the "drought subset" (any one
   of the K drought objectives above its historical median is enough to
   call a trace droughty) and the **intersection complement** for the
   "nominal subset" (all K objectives below the historical median).
   Justification: the union generously admits traces that achieve
   drought intensity on at least one axis, which is the relevant test
   for the FDC-low-flow criterion; the strict intersection complement
   keeps the nominal subset hydrologically unremarkable, which is the
   relevant test for the "does the archive still contain historical-
   looking traces" criterion. Revise if the resulting subsets have
   fewer than 100 solutions at any relevant NFE.
3. Should Criterion 4's "70 % of historical spread" threshold be the
   bar, or do we want the synthetic archive to match the historical
   spread exactly, or exceed it? **A: Match the spread exactly.**
   Criterion 4 pass test updated to ``0.9 ≤ std(A_syn) / std(A_hist) ≤ 1.1``.
4. What is the minimum number of Pareto solutions per "subset" for
   Criteria 2-3 to be meaningful? **A: Target > 100 solutions per
   subset.** If either the drought_subset or nominal_subset has < 100
   solutions, the envelope comparison is flagged unreliable and the
   criterion is reported as "insufficient data, ramp NFE" rather than
   pass/fail.
