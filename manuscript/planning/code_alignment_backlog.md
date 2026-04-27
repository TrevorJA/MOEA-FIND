# Code Alignment Backlog

*Originally prepared 2026-04-15 from comparison of `manuscript_main_draft.md` §2.2.4 against `src/constraints.py` and the analytic experiment scripts. Bootstrap-calibration guidance referenced below as `generator_constraints_review.md` was rolled into per-experiment SI reporting and the diagnostic at `workflows/diagnostics/diag_constraint_calibration.py`.*

*Framing principle: the manuscript is the specification. Items below are places where the code has not yet been updated to match it. Items marked [TREV-DECISION] require authorial choice before the code alignment can be completed.*

*2026-04-19 update: Items 1, 2, and 4 are retired. The manuscript §2.2.4 no longer pre-specifies constraint tolerances. Instead, tolerances are calibrated by bootstrap per experiment (via `workflows/diagnostics/diag_constraint_calibration.py`) and reported in the Supporting Information for each run. The pre-2026-04-19 §2.2.4 prose hard-coded lag-1 AC = 0.05 and non-drought mean = 15 %, which the 200K NFE Cannonsville run did not use; softening the prose to the per-experiment-reporting posture eliminates the mismatch without forcing a code-side change. Items 3, 5, 6, 7 remain active.*

---

## Item 1: Lag-1 autocorrelation tolerance

**Manuscript** (§2.2.4, paragraph 2):
> "The lag-1 autocorrelation of the trace is required to lie within $0.05$ of the historical lag-1 autocorrelation..."

**Current code** (`src/constraints.py`, `autocorrelation_constraint`, default `tolerance=0.3`; `compute_all_constraints`, default `autocorr_tol=0.3`):
The effective tolerance at runtime is 0.30. No call site overrides this default.

**Context:** The constraint verifies statistical plausibility — that the synthetic trace exhibits lag-1 autocorrelation consistent with what the Kirsch generator could produce given an infinite ensemble. The manuscript value of 0.05 reflects a tighter plausibility envelope than the current code. `generator_constraints_review.md` §4.4 estimates a bootstrap-calibrated plausibility tolerance of 0.08–0.12 for 30-year traces.

**Alignment action:** Update `autocorr_tol` default in `compute_all_constraints` to the chosen value and verify no call site requires the wide 0.30 tolerance. [TREV-DECISION: final tolerance value — manuscript states 0.05; calibrated range per §4.4 is 0.08–0.12.]

---

## Item 2: Non-drought period mean constraint not implemented

**Manuscript** (§2.2.4, paragraph 2):
> "...and the mean of the non-drought annual flows is required to lie within $15$ percent of the historical average."

**Current code** (`src/constraints.py`, `annual_stats_constraint`):
```python
s_mean = np.mean(synthetic_monthly)   # mean over ALL months
h_mean = np.mean(historical_monthly)
mean_violation = max(0.0, abs(mean_ratio - 1.0) - mean_tolerance)  # tolerance = 0.50
```
The function computes the mean over every monthly value, not over non-drought months only. The tolerance is 0.50, not 0.15. No drought-period filtering exists in `constraints.py`.

**Context:** The manuscript describes a plausibility check on non-drought periods specifically — verifying that periods between drought events remain statistically consistent with the generator's historical envelope. This is a narrower and more interpretable plausibility criterion than the all-flow version. Implementing it requires per-trace SSI event detection to identify non-drought months.

**Alignment action:** Implement the non-drought-period mean constraint as described in the manuscript, following `generator_constraints_review.md` §4.3. [TREV-DECISION: confirm final tolerance value — manuscript states 15%; calibrated guidance is in §4.3.]

---

## Item 3: Seasonal cycle constraint not described in manuscript

**Manuscript** (§2.2.4, paragraph 2):
> "Two plausibility constraints are imposed on generated traces. The lag-1 autocorrelation...and the mean of the non-drought annual flows..."

**Current code** (`src/constraints.py`, `seasonal_cycle_constraint` and `compute_all_constraints`):
A third constraint, `seasonal_cycle_constraint` (default `tolerance=0.5`), checks that per-month means of the synthetic trace lie within 50 percent of the historical per-month means for all 12 calendar months.

**Context:** The seasonal cycle constraint at 50 percent tolerance is approximately redundant given that the Kirsch generator preserves monthly geometric means exactly by construction (§2.2.1). Its presence does not change results but its absence from the manuscript description is an inaccuracy.

**Alignment action:** Add one sentence to §2.2.4 naming the seasonal cycle constraint and noting its near-redundancy with the generator's exact mean preservation. No TREV-DECISION required.

---

## Item 4: Annual mean and CV tolerance — undescribed constraints

**Manuscript** (§2.2.4): Neither the annual mean constraint nor the CV constraint is mentioned.

**Current code** (`src/constraints.py`, `annual_stats_constraint`, defaults `mean_tolerance=0.5`, `cv_tolerance=0.5`):
Both annual mean and annual CV of synthetic monthly flows must lie within 50 percent of their historical counterparts.

**Context:** `generator_constraints_review.md` §4.1 flags the 50 percent tolerance as not defensible and recommends a bootstrap-calibrated tolerance near 0.15. The manuscript should describe these constraints, and the code should use the chosen tolerance.

**Alignment action:** Add the annual mean and CV constraints to the §2.2.4 description. Update code tolerances to the chosen values. [TREV-DECISION: final tolerance values — current code uses 0.50; calibrated recommendation per §4.1 is 0.15.]

---

## Item 5: Analytic benchmark scripts use EpsNSGAII — note only, not a discrepancy

**Manuscript** (§2.2.2):
> "The multi-objective evolutionary optimizer used by MOEA-FIND is Borg MOEA (Hadka and Reed, 2013)..."

**Current analytic scripts** (`scripts/02_analytic_3d.py`, `scripts/01_analytic_2d.py`):
```python
from platypus import EpsNSGAII, Problem, Real
algo = EpsNSGAII(problem, epsilons=[epsilon] * (K + 1))
```
The analytic benchmark scripts use EpsNSGAII from `platypus` as a locally runnable stand-in. The single-site Cannonsville script (`04_kirsch_single_site.py`) docstring reads "serial Borg/platypus fallback," confirming that platypus is the local substitute and Borg is the production algorithm.

**This is not a manuscript discrepancy.** Borg is the production algorithm; EpsNSGAII is the local test stand-in used for rapid development on machines without the Borg license. The manuscript correctly identifies Borg. The analytic results are valid because both algorithms use the same epsilon-dominance archive mechanism, which is the only property the theoretical argument requires. The Section 3.1 dimension sweep figures produced with EpsNSGAII are valid preliminary results pending reproduction with Borg on HPC.

**Note for §3.1:** A parenthetical identifying EpsNSGAII as the algorithm used for the analytic runs, with a note that production Cannonsville runs use Borg, cleanly separates preliminary from production results without contradicting the §2.2.2 algorithm description.

---

## Item 6: NFE for analytic benchmark — two runs conflated

**Manuscript** (§3.1): "with 30,000 function evaluations per $K$"

**Data files:**
- `outputs/diag_shell_vs_interior/k3/config.json`: `nfe=30000` (Figure 4 dimension sweep data)
- `outputs/exp02_analytic_3d/config.json`: `nfe=2000` (Supporting Information SI-1 hyperplane figure data)

**Context:** The manuscript §3.1 correctly states 30,000 NFE for Figure 4. However, the affine identity verification in SI-1 uses the 2,000-NFE run (682 archive members), not the 30,000-NFE, 6,158-member Figure 4 run.

**Alignment action:** Add a parenthetical in §3.1 clarifying that the affine identity verification in SI-1 uses the 2,000-NFE run (682 archive members) and the coverage diagnostics in Figure 4 use the 30,000-NFE dimension sweep. No TREV-DECISION required.

---

## Item 7: CLAUDE.md archive size — outdated session note

**CLAUDE.md** (project notes): "1362 Pareto solutions, hyperplane verified to machine precision (std~10^-16)"

**Data files:**
- `outputs/diag_shell_vs_interior/k3/results.json`: n=6,158 at K=3
- `outputs/exp02_analytic_3d/results.json`: n=682 at K=3

**Context:** The 1,362 figure appears to derive from an intermediate exploratory run no longer present in the outputs directory. The revised §3.1 uses the actual output-file numbers (6,158 for Figure 4; 682 for SI-1).

**Alignment action:** Update CLAUDE.md to reflect the authoritative output-file numbers from the current runs. Not a manuscript alignment item but worth tracking here to prevent the stale figure from re-entering the manuscript.

---

## Summary table

| # | Location | Manuscript specification | Current code | Action |
|---|---|---|---|---|
| 1 | §2.2.4 | lag-1 AC plausibility tol = 0.05 | tol = 0.30 | Update code; TREV-DECISION on final value |
| 2 | §2.2.4 | non-drought mean within 15% | all-flow mean within 50% | Implement non-drought filter in code; TREV-DECISION on final value |
| 3 | §2.2.4 | two constraints listed | three present (seasonal cycle also) | Add one sentence to manuscript describing seasonal cycle constraint |
| 4 | §2.2.4 | annual mean/CV not described | annual mean + CV at 50% tol | Describe in manuscript; update code tolerances; TREV-DECISION on final values |
| 5 | §2.2.2, §3.1 | Borg MOEA (production) | EpsNSGAII (local stand-in) | No change needed — note in §3.1 parenthetical only |
| 6 | §3.1 | implies single NFE for all analytic results | SI-1 from 2K run; Fig 4 from 30K run | Add parenthetical clarifying the two runs |
| 7 | CLAUDE.md | — | 1362 (stale) vs 6158/682 (current outputs) | Update CLAUDE.md |

---

## Item 8: Tighten anti-ideal `D*` in Script 04 optimization runs

**Issue (2026-04-20):** The anti-ideal point read from
`outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json`
(`"anti_ideal"`) sits far outside the realized drought-characteristics
cloud from the 500-member Stage 1 replay. When plotted on fig09 it
dominated the axes and pushed the data into a corner; fig09 now crops
the anti-ideal off-plot, but the underlying issue is that `D*` defines
the Manhattan target the MOEA optimizes toward, so an overly loose `D*`
produces a correspondingly loose feasible region.

**Alignment action:** In the next `04_kirsch_single_site.py`
optimization run, reduce `D*` toward the empirical worst-case
duration/severity observed in Stage 1 generators (or a calibrated
percentile such as the 95th of the generated ensemble). This tightens
the Manhattan objective and should bring the anti-ideal inside the
fig09 axes without needing to clip it. [TREV-DECISION: final `D*`
value — empirical worst-case vs fixed percentile.]

---

## Item 9: NYC demand as a scenario-discovery dimension

**Idea (2026-04-20):** Factorially sample NYC/NJ demand across a range
(e.g., 400–800 MGD, or a multiplicative scaling factor around the
FFMP baseline) and pair every demand level with every Pareto drought
scenario. Scenario discovery then operates in 4D feature space
(duration, severity, peak_month, demand), letting the classifier
learn demand effects separately from hydrology.

**Precedent:** Hadjimichael et al. 2020, Kasprzyk & Reed (and related
Colorado River SD work) use multiplicative demand factors as
uncertainty dimensions when doing SD for water-resources systems.
Treating demand as a "deeply uncertain" axis alongside hydrology is
the convention in that literature.

**Why:** Under the current configuration, demand is held at the FFMP
`max_flow_baseline_delivery_nyc` scalar — an extreme-demand
stress-test. Real demand uncertainty over a policy horizon spans a
meaningfully wider range; if the classifier learns "all observed
non-satisficing scenarios have severity > 1.5", that boundary could
shift significantly under lower demand. Having demand as an explicit
dimension makes those sensitivities visible.

**How to apply:** Defer until after the core MOEA-FIND paper lands,
then run as a separate numerical experiment. Implementation sketch:

- Add `--demand-scales 0.7 0.85 1.0 1.15 1.3` (for example) to
  `09_drb_policy_reeval.py`. For each scale, run
  `N_demand × N_pareto` Pywr-DRB simulations. Either patch Pywr-DRB
  ModelBuilder to accept an explicit scalar override (instead of the
  hardcoded `max_flow_baseline_delivery_nyc` lookup) or generate a
  `custom` flow_type per scale.
- Extend the metric bank schema to carry the per-scenario demand scale
  as a column.
- Fig09 extension options: (a) one row of panels per demand scale —
  panels still 3D-sliced on drought characteristics, demand held at
  the row's scale; (b) single-row panels but color-code scatter by
  demand scale; (c) a separate supplementary figure stratified by
  demand.

**Challenges:**

1. **Compute cost:** Simulation cost multiplies by `N_demand`. Five
   scales × 500 Pareto members = 2500 pywrdrb runs. At ~75 min per
   batch of 500 that's ~6 h of Stage 3 time — manageable but non-trivial.
2. **Plotting in 4D:** fig09 already stretches a 3D feature space
   onto 2D panels; demand as a 4th axis forces either more panels or
   a faceting grid. Risk of clutter.
3. **Narrative scope:** If the MOEA-FIND manuscript is primarily a
   *methodological* contribution about synthetic drought generation,
   demand-dimension SD may dilute the story. Better suited to a
   follow-up paper that uses MOEA-FIND's outputs as one input among
   several to a deeply-uncertain-futures analysis.

**Alignment action:** No code change yet. Revisit post-submission of
the current manuscript. [TREV-DECISION: paper scope — include demand
dimension in the core paper's fig09 or relegate to a follow-up?]

---

## Item 10: Update CLAUDE.md to reflect post-2026-04-20 HPC pull

**Issue (2026-04-27):** `.claude/CLAUDE.md` is dated 2026-04-13 (Session
14). The 2026-04-20 and 2026-04-21 HPC pulls added seven production
modules to `src/` (`pywrdrb_bridge.py`, `multisite_data.py`,
`scenario_discovery.py`, `satisficing_metrics.py`, `satisficing_labels.py`,
`constraints_dv.py`, `historical_blocks.py`) and reorganised experiments
into `workflows/experiments/` (numbered 01–17) with one-to-one SLURM
wrappers under `workflows/slurm/`. None of this is reflected in
CLAUDE.md, including the explicit "Pending Phase 2 Experiment 2.1, 2.3"
note which is no longer accurate.

**Alignment action:** Replace the "Current State" section and the
"Project Structure" tree in `.claude/CLAUDE.md` with the snapshot from
`code_state.md`. Update the immediate-next-steps list. Defer until after
the next HPC results land so the CLAUDE.md update is a single sweep
rather than two.

---

## Item 11: Manuscript notes on the HPC integration modules

**Issue (2026-04-27):** The seven HPC integration modules referenced in
Item 10 are documented only in their docstrings. The Methods section of
the manuscript currently describes Sections 3.2 and 3.3 in terms of
"Pywr-DRB simulation" without naming the bridge architecture. The Methods
prose may need a paragraph or sub-subsection covering the three-stage
bridge (replay → preprocess → simulate) and the satisficing-label fixed
panel.

**Alignment action:** When Section 3.3 prose is drafted, add a paragraph
in §2 or §3 describing the bridge stages and the satisficing metric
panel, with citations to the FFMP drought-level documentation and the
Hashimoto et al. (1982) reliability metrics. No code change. [TREV-DECISION:
whether the bridge architecture belongs in Methods or Supporting
Information.]
