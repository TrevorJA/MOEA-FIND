# Statistical Constraint Review for MOEA-FIND

*Prepared 2026-04-14. Addresses reviewer critique 2 (Section 2.2.4 of manuscript draft).
Covers: (1) Kirsch-Nowak preservation analysis from generator source, (2) taxonomy
of standard validation statistics in stochastic hydrology, (3) gap analysis of current
MOEA-FIND constraints, (4) recommended constraint set with calibration rationale,
(5) proposed SI sensitivity analysis, (6) priority ranking.*

---

## 1. Background: What the Kirsch-Nowak Generator Preserves by Construction

The Kirsch-Nowak generator (Kirsch et al., 2013; Nowak et al., 2010) as implemented
in SynHydro (`kirsch.py`, v2026) follows a seven-step pipeline. Understanding which
statistics each step fixes, approximates, or ignores is prerequisite to any constraint
design, because adding a constraint that enforces a statistic the generator already
fixes exactly is computationally wasteful and contributes no plausibility protection.

### 1.1 The seven-step pipeline

**Step 1 (log transform and standardization).** Each observed monthly flow is
log-transformed and standardized by the calendar-month mean and standard deviation
in log-space. This produces the standardized residual matrix Z_h of shape (n_years,
12, n_sites), stored in `self.Z_h`. The per-month log-space mean vector `self.mean_month`
and standard deviation vector `self.std_month` are computed from the historical record
and stored.

**Step 2 (normal-score transform).** When `generate_using_log_flow=True` (the default),
the residuals Z_h are mapped through the Hazen plotting-position empirical CDF for each
(month, site) combination to produce approximately normal scores Y. This step maps each
calendar-month marginal distribution in log-standardized space to N(0,1) regardless of
the shape of the empirical distribution. The inverse transform (Step 6) maps back through
the same empirical CDF with linear extrapolation at the tails.

**Step 3 (shifted Y_prime for Dec-Jan boundary).** A second matrix Y_prime is formed
by shifting Y forward by six months so that the first six months of Y_prime[i] correspond
to July--December of year i and the last six months correspond to January--June of year
i+1. This construction allows the Cholesky decomposition (Step 4) to impose correlation
across the December-to-January calendar boundary.

**Step 4 (Cholesky decomposition and application).** Two 12-by-12 correlation matrices
are computed: one from Y over all years (stored as `U_site[s]`) and one from Y_prime
(stored as `U_prime_site[s]`). Both are decomposed via Cholesky factorization. Bootstrap
draws are decorrelated by simple indexing and then correlated by matrix multiplication
with the Cholesky factor: `Z[:,:,s] = X[:,:,s] @ U_site[s]`. This operation imposes
the historical 12-by-12 monthly correlation structure in normal-score space.

**Step 5 (combine Z and Z_prime).** The first six months of the synthetic year come
from the Z_prime branch (carrying Dec-Jan cross-year information) and the second six
months come from the Z branch. This combination is the `_combine_Z_and_Z_prime` operation.

**Step 6 (inverse normal-score transform).** The combined Z is mapped back through
the inverse empirical CDF for each (month, site). The tails are linearly extrapolated
so that synthetic values drawn from outside the historical quantile range are not clamped
to the historical min/max.

**Step 7 (destandardization and exponentiation).** Synthetic flows are recovered by
reversing Step 1: multiply by `std_month`, add `mean_month`, then exponentiate.

### 1.2 Statistics preserved exactly by construction

The following statistics are algebraically fixed by the generator design regardless
of which bootstrap indices or residuals are injected.

**Monthly geometric means.** Step 7 adds back `mean_month` (the historical log-space
per-month mean) and multiplies by `std_month`, then exponentiates. For any bootstrap
draw that passes through this pipeline, the expected value of the synthetic log-flow
for calendar month m equals the historical log-space mean for month m. Over a long
trace, the sample mean of the log-transformed flows for each calendar month converges
to the historical value. This makes the monthly geometric mean an exact asymptotic
invariant of the generator. For a 30-year trace, finite-sample variance around this
mean is on the order of (historical log-space std)/sqrt(30).

**Monthly log-space standard deviation.** Analogous argument: the rescaling by
`std_month` in Step 7 imposes the historical log-space standard deviation per
calendar month by construction.

**Within-year monthly correlation structure.** The Cholesky multiplication of Step 4
imposes the historical 12-by-12 correlation matrix exactly in normal-score space.
The pairwise correlation between any two calendar months within a single water year,
computed from a large synthetic ensemble, will converge to the historical value in
normal-score space. Back-transformation through the empirical CDF and exponentiation
introduces a monotone transform, so the rank correlation between any two within-year
months is also preserved exactly; the linear correlation coefficient is approximately
preserved to the extent that the joint distribution is close to bivariate normal in
normal-score space.

**December-to-January lag correlation (approximately).** The Y_prime construction
preserves the Dec-Jan cross-year correlation in the same manner as the within-year
structure. This is an approximation rather than an exact fix because only the boundary
six-month blocks are handled; correlations at lags of two or more calendar years are
not explicitly imposed.

**Flow non-negativity.** The exponentiation in Step 7 guarantees Q > 0 for all
synthetic values. This cannot be violated regardless of the decision vector injected
by the MOEA.

**Per-month empirical marginal distribution.** The combination of normal-score transform
(Step 2) and inverse transform (Step 6) constitutes a nonparametric probability integral
transform. When a bootstrap index is used for a given calendar month, the residual drawn
is one of the n_years historical residuals for that month, which is exactly rank-order
consistent with the historical marginal. For residual-mode injection, the DV p is
mapped through `self.sorted_residuals[idx, month, site]` where idx is the floor-quantile
of p, so the residual is drawn from the historical empirical CDF. The monthly flow
duration curve (FDC) per calendar month is therefore approximately preserved by the
generator's sampling machinery; it is not exactly preserved at finite trace length
but converges to the historical FDC as trace length increases.

### 1.3 Statistics the generator does not preserve by construction

The following statistics are not fixed by any step in the pipeline. These are the
degrees of freedom that MOEA-FIND exploits, and any of them can be pushed away from
historical values by the optimizer without triggering a violation of the generator's
internal logic.

**Inter-annual autocorrelation and Hurst exponent.** The bootstrap in the standard
generator draws year indices independently at each calendar month position. Year-to-year
dependence is introduced only at the Dec-Jan boundary through Y_prime; there is no
constraint on correlations at lags of 2, 3, or more years. Streamflow records in
humid basins typically exhibit Hurst exponents H in the range 0.55 to 0.80
(Montanari et al., 1997; Salas, 1993), reflecting slow decorrelation of multi-year
storage states. The i.i.d. bootstrap gives H = 0.50 asymptotically. When MOEA-FIND
selects a sequence of historically low-flow year indices to construct a prolonged
drought, it implicitly creates positive inter-annual autocorrelation that may inflate
H beyond the historical value. Conversely, selection of alternating low-flow and
high-flow years to meet an intermediate drought objective can suppress apparent
persistence below the historical value. Wheeler et al. (2024) explicitly include the
Hurst coefficient as an optimization objective in their multisite extension of the
Borgomeo search-based generation framework for precisely this reason.

**Annual flow volume and annual CV.** The annual total flow of a synthetic trace is
the sum of its twelve monthly flows and is not fixed by the monthly-mean constraints.
A trace in which the optimizer selects 30 consecutive historically low-flow year
indices can have an annual mean well below the historical value (potentially 50% below
for extreme drought traces), while still satisfying the monthly correlation and
seasonal cycle constraints. The current code imposes `mean_tolerance=0.5` and
`cv_tolerance=0.5` in `annual_stats_constraint`, which allows a factor-of-two
excursion in either direction around the historical annual mean. This is an extremely
wide band by any reasonable hydrological standard.

**High and low flow quantiles (annual FDC).** The distribution of annual total flows
across the 30-year trace is not constrained beyond the annual mean and CV (which are
themselves loosely constrained). The 10th and 90th percentile of the annual flow
distribution can deviate substantially from historical without violating any current
constraint, because the current constraint operates on the trace mean rather than on
quantile structure.

**Lag-1 autocorrelation of the full monthly series (effective constraint vs. stated
tolerance).** The manuscript (Section 2.2.4) states a tolerance of 0.05 for the lag-1
autocorrelation constraint, but the current implementation in `constraints.py` uses
`autocorr_tol=0.3`. These differ by a factor of six. At the actual implemented tolerance
of 0.30, the constraint is effectively inactive for the historical Cannonsville record,
whose lag-1 autocorrelation is approximately 0.85 (a value close to 1.0 is typical
for monthly streamflow due to seasonal cycle). A deviation of 0.30 from 0.85 allows
lag-1 values as low as 0.55, which would represent a physically unusual streamflow
record. The stated 0.05 tolerance would be appropriate but is not what is implemented.
This discrepancy must be resolved before submission.

**Monthly skewness.** The log-normal assumption embedded in Step 1 implicitly constrains
skewness in log-space to zero, but the skewness of flows in physical space depends on
the magnitude of the log-space variance. For a targeted trace in which the optimizer
selects residuals consistently from the lower tail of the monthly distribution, the
skewness of physical-space flows can shift from the historical value. No constraint
currently addresses this.

**Non-drought period flows (implementation gap).** The manuscript states that "the mean
of the non-drought annual flows is required to lie within 15 percent of the historical
average." The actual implementation (`annual_stats_constraint` in `constraints.py`)
computes the mean over all monthly flows, not over non-drought months. The drought-
period/non-drought-period split is event-detected per trace and is computed only in the
objectives pipeline, not in the constraint. The constraint as implemented does not match
the manuscript claim and does not prevent the optimizer from depressing flows across the
full record rather than only during drought periods.

---

## 2. Taxonomy of Standard Validation Statistics in Synthetic Hydrology

The canonical validation checklist for stochastic streamflow generators is drawn from
Stedinger (1982), Salas et al. (1982), and the synthetic scenario generation literature
initiated by Borgomeo et al. (2015) and extended by Zaniolo et al. (2023) and Wheeler
et al. (2024). Authors presenting a new generator are expected to demonstrate reproduction
of statistics in each of the following five categories. Deficit in any category is
grounds for reviewer concern.

### 2.1 Marginal statistics

- **Monthly mean flow** (arithmetic or geometric, per calendar month): verified by
  visual comparison of box plots or violin plots across synthetic ensemble vs. historical.
- **Monthly standard deviation** (or CV) per calendar month.
- **Monthly skewness coefficient** per calendar month: particularly important when the
  generator uses a log-normal assumption, because real streamflow skewness departs from
  log-normal in low-flow months.
- **Annual mean flow** and **annual CV**: aggregate statistics that summarize the
  trace-level water balance.
- **Extreme quantiles** of the annual distribution: P5, P10, P90, P95 of annual totals
  or of monthly flows pooled over all calendar months. These directly affect drought
  and flood characterization.

### 2.2 Temporal dependence

- **Lag-1 autocorrelation** of monthly flows: the foundational check in ARMA-based
  hydrology (Thomas and Fiering, 1962) and retained in all modern generators.
- **Lag-2 autocorrelation**: sensitive to multi-month drought clusters that lag-1 alone
  does not capture.
- **Hurst exponent** (H): the single most important summary of long-range temporal
  dependence in streamflow. H = 0.5 corresponds to i.i.d. structure (white noise),
  H > 0.5 to positive long-memory (clustered low-flow and high-flow periods), and H < 0.5
  to anti-persistent structure. Streamflow in humid watersheds typically has H between
  0.55 and 0.80 (Montanari et al., 1997). Violation of H matching is cited as a
  known failure mode of purely i.i.d. bootstrap generators (Salas, 1993). Wheeler et al.
  (2024) and Borgomeo et al. (2015) both flag Hurst preservation as a key validation
  criterion for search-based generators.
- **Spectral density** at low frequencies: related to the Hurst exponent but expressed
  in the frequency domain. Not typically required for a first submission but may be
  requested by reviewers with a long-memory background.

### 2.3 Cross-site spatial dependence

- **Cross-site Pearson correlation matrix** of concurrent monthly flows: the most common
  spatial check. Kirsch preserves this by construction for the index-injection mode
  (shared index across sites); it must be verified empirically for residual-injection
  mode because the per-DV mapping does not guarantee identical rank selection across sites.
- **Cross-site Spearman rank correlation**: more robust to marginal distribution shape,
  preferred when monthly distributions are highly skewed.
- **Upper and lower tail dependence between sites**: relevant when the analysis targets
  simultaneous multi-site droughts. Copula-based generators (Brechmann et al., 2017)
  are specifically designed to preserve this; bootstrap generators only approximately
  do so through the shared index.

For the single-site Cannonsville application, cross-site statistics are not directly
relevant and can be acknowledged as out-of-scope for the present paper.

### 2.4 Drought and extreme statistics

- **Run-length distribution** of below-threshold events: distribution of drought
  durations, not just mean duration.
- **Deficit volume distribution**: distribution of cumulative deficits per event.
- **Return period of multi-year droughts**: for a 30-year trace, the 10-year and 20-year
  return period droughts should remain physically plausible.
- **Maximum monthly flow and minimum monthly flow**: bounds check, important for
  confirming no unrealistic extremes are generated.
- **SSI distribution**: for SSI-based drought objectives, the distribution of SSI values
  from the generator should match the reference distribution (Gamma or log-normal fitted
  to historical) across the full range, not just at the drought threshold.

### 2.5 Mass balance and hydrologic plausibility

- **Non-negative flows**: guaranteed by the generator.
- **Seasonal cycle (monthly mean pattern)**: verified by the existing seasonal cycle
  constraint in MOEA-FIND.
- **Annual water balance (annual total flows within historical range)**: partially covered
  by the annual stats constraint but at an excessively wide tolerance.
- **Absence of physically impossible patterns**: flows that are simultaneously zero and
  positive, or that exhibit negative lag-1 autocorrelation in a basin with strong
  seasonal storage. These are typically filtered by the generator design but can arise
  from adversarial index selection.

---

## 3. Gap Analysis: Current MOEA-FIND Constraints

The table below records, for each major statistic in Section 2, whether the Kirsch-Nowak
generator fixes it by construction, whether the current MOEA-FIND constraint set covers
it, and whether the coverage is adequate for submission.

| Statistic | Fixed by Kirsch | Covered by MOEA-FIND | Adequate? | Notes |
|---|---|---|---|---|
| Monthly geometric mean | Yes (exact) | Redundant | Yes | No constraint needed |
| Monthly log-space std | Yes (exact) | Redundant | Yes | No constraint needed |
| Within-year monthly correlations | Yes (exact) | Redundant | Yes | No constraint needed |
| Dec-Jan lag correlation | Yes (approx) | Redundant | Yes | No constraint needed |
| Flow non-negativity | Yes (exact) | Redundant | Yes | No constraint needed |
| Per-month empirical FDC | Yes (approx) | Redundant | Yes | Approximately redundant |
| Lag-1 AC of monthly series | No | Yes (tol=0.30) | No | Manuscript says 0.05; code says 0.30 |
| Lag-2 AC of monthly series | No | No | No | Not covered |
| Hurst exponent | No | No | No | Critical gap; Wheeler 2024 includes this |
| Annual mean flow | No | Yes (tol=0.50) | No | 50% tolerance is far too wide |
| Annual CV | No | Yes (tol=0.50) | No | Same |
| Annual flow quantiles (P10, P90) | No | No | No | Not covered |
| Monthly skewness | No | No | Marginal | Low priority; exp() partially constrains |
| Non-drought period flows | No | No (code mismatch) | No | Manuscript claim not implemented in code |
| Seasonal cycle (monthly means) | No | Yes (tol=0.50) | No | Tolerance calibration needed |
| Cross-site correlation | Yes for index mode | Not explicitly | Partial | Multi-site extension only |

### 3.1 Specific unconstrained degrees of freedom the MOEA can exploit

The optimizer can produce a "mathematically optimal but hydrologically implausible"
trace through three distinct mechanisms that the current constraint set does not prevent.

**Mechanism 1: Temporal clustering of identical or near-identical historical years.**
In index-injection mode, nothing prevents the optimizer from selecting the same drought
year index repeatedly. A trace composed of 30 copies of the worst historical drought
year would have a statistically valid within-year correlation structure (because the
Cholesky factor is applied at each year separately) and would pass the seasonal cycle
constraint (because the monthly means reproduce the historical drought year's monthly
pattern, which is some fraction of the historical long-term monthly means). Yet such
a trace is physically implausible and not a plausible representation of a future 30-year
period. The lag-1 autocorrelation constraint catches some but not all manifestations
of this, because repeatedly selecting the same year produces a nearly constant series,
which has lag-1 AC close to 1.0 rather than the historical value near 0.85 for monthly
flows. At the implemented tolerance of 0.30, this would be flagged; at any wider
tolerance, it would not.

**Mechanism 2: Depressing the full-trace annual mean under the guise of drought
statistics.** The optimizer's drought objectives reward low SSI values, which correlate
with low annual flows. A trace that depresses all monthly flows by 40% below historical
would have extreme drought characteristics but would still pass the annual stats constraint
at the current 50% tolerance. The 15% non-drought period constraint stated in the
manuscript is not implemented, so this mechanism is entirely unconstrained in the code.

**Mechanism 3: Inflating apparent Hurst exponent through adversarial year-sequence
selection.** The optimizer can select sequences of low-flow years for the drought period
and high-flow years for the non-drought period, creating a biphasic trace with higher
inter-annual autocorrelation than the historical record. This inflates H above the
historical value. Conversely, selecting alternating low-flow and high-flow years at
annual scale creates anti-persistent structure inconsistent with streamflow physics.
Neither behavior is caught by any current constraint.

---

## 4. Recommended Constraint Set

For each weak point identified in Section 3, I propose a specific constraint with
calibration rationale, computational cost, and implementation recommendation.

### 4.1 Annual water balance constraint (replace current annual stats constraint)

**Statistic constrained.** Mean annual flow volume (arithmetic mean of annual flow
totals across the trace), reported as a fractional deviation from the historical mean.

**Current implementation.** `annual_stats_constraint` in `constraints.py` with
`mean_tolerance=0.5` and `cv_tolerance=0.5`, applied to the full monthly flow vector
without drought-period filtering.

**Proposed replacement.** Constrain the mean annual flow to lie within the 95%
bootstrap confidence interval of the historical mean annual flow, computed by
resampling the historical record with replacement at annual resolution.

For the Cannonsville record (water years 1952-2022, N=70 years), drawing 30-year
subsamples with replacement and computing the sample mean annual flow for each draw
produces a 95% CI approximately equal to the historical mean plus or minus one
standard error times 1.96, where the standard error is the historical annual flow
standard deviation divided by sqrt(30). For a typical CV of 0.30 on annual flows,
this gives a 95% CI width of approximately plus or minus 0.11 of the historical mean.
A symmetric tolerance of 0.15 is therefore slightly wider than the 95% bootstrap CI
and is the appropriate value for this constraint.

**Tolerance.** 0.15 (fractional deviation from historical mean annual flow). This
replaces the current 0.50 tolerance. The 0.15 value is calibrated to the 95% bootstrap
CI of the 30-year sample mean given historical annual flow variability; it is not a
round number chosen arbitrarily.

**Computational cost.** One sum and one division per trace evaluation: O(n_months).
Negligible.

**Constraint type.** Hard constraint (infeasibility flag) for deviations beyond 0.30
(twice the target tolerance); soft constraint (penalty proportional to excess deviation)
between 0.15 and 0.30. This two-tier implementation prevents the optimizer from
concentrating on traces near the constraint boundary while allowing exploration of the
full feasible drought space.

**Additional note.** The CV constraint should be retained at a similarly calibrated
tolerance. Historical annual flow CV is approximately 0.30 for Cannonsville. Allowing
the synthetic trace CV to deviate by more than 0.30 from the historical CV (i.e.,
current 50% tolerance on the CV ratio) permits traces whose inter-annual variability
is double or half the historical value, which represents a qualitatively different
hydroclimatic regime. A tolerance of 0.30 on the CV ratio (fractional deviation) is
recommended, corresponding to the approximate 95% bootstrap CI of the sample CV at
30-year trace length given a population CV of 0.30.

### 4.2 Hurst exponent constraint

**Statistic constrained.** The Hurst exponent H of the monthly flow series, estimated
by the rescaled range (R/S) analysis (Hurst, 1951) or detrended fluctuation analysis
(DFA; Peng et al., 1994). DFA is preferred because it is less sensitive to the choice
of record length subdivision than R/S and is directly applicable to monthly flow series
of 360 observations (30 years times 12 months).

**Rationale.** The Kirsch i.i.d. bootstrap gives H approaching 0.50 for large
ensembles; the target historical value for Cannonsville should be computed empirically
from the historical record. Wheeler et al. (2024) include Hurst coefficient as an
explicit optimization objective precisely because their simulated annealing search
distorts it. MOEA-FIND faces the same distortion risk when the optimizer selects
adversarial year sequences to maximize drought characteristics.

**Tolerance calibration.** Bootstrap the historical Cannonsville record: draw 1000
subsets of 30 years with replacement (to match synthetic trace length) and compute
the DFA-estimated H for each. The resulting empirical distribution gives the 95%
CI for H at trace length T=360 months. For streamflow records of typical persistence,
this CI is approximately plus or minus 0.15 around the historical H. The constraint
tolerance should be set to the half-width of this bootstrap CI.

**Computational cost.** DFA requires O(n_months * log(n_months)) operations per trace
evaluation, approximately 360 * 9 = 3240 floating-point operations. At typical function
evaluation rates for a Kirsch trace (order 10 milliseconds per evaluation), the DFA
computation adds less than 1% overhead. This is acceptable.

**Constraint type.** Soft constraint (penalty proportional to excess deviation outside
the bootstrap CI bounds). Hurst exponent estimation at 360 observations has substantial
finite-sample variance, making a hard infeasibility boundary inappropriate. A quadratic
penalty that activates when H_synthetic falls outside the 95% bootstrap CI is
recommended.

**Implementation note.** The DFA function must be applied to the full monthly flow
sequence, not to the annual flow sequence (which has only 30 observations, too few for
reliable Hurst estimation). Log-transformation before DFA is optional but recommended
for right-skewed series.

### 4.3 Non-drought period mean flow constraint (implement as stated in manuscript)

**Statistic constrained.** Mean flow across all months classified as non-drought months
in the SSI-3 event extraction applied to the synthetic trace.

**Current implementation gap.** The manuscript states this constraint at 15%
tolerance applied to non-drought period flows, but the code implements it as a
constraint on the mean of all flows (including drought months) at 50% tolerance.
These are not the same constraint. The 15% non-drought figure must be implemented
as stated or the manuscript claim must be corrected.

**Clarification needed: fixed vs. event-detected split.** The reviewer critique
(critique 2 in reviewer_critiques.md) explicitly asks whether the split is computed
on a fixed calendar partition or on the event-detected split of each trace. The answer
must be event-detected per trace, because the drought period length changes as the
optimizer generates more extreme traces. A fixed calendar split would make the
constraint blind to traces whose drought occupies unusual months.

**Tolerance calibration.** The 15% figure should be calibrated against the bootstrap
distribution of non-drought period mean flows computed from the Kirsch library.
Generate 1000 standard Kirsch traces without optimization, extract the non-drought
months using the same SSI-3 event detection algorithm, and compute the mean non-drought
flow for each trace. The 95% CI of this distribution gives the appropriate tolerance.
For a generator that approximately preserves the historical FDC, the non-drought flow
mean should be close to the historical mean; its 95% CI at 30-year trace length is
roughly plus or minus 10-15% of the historical value. The 15% figure is therefore
approximately appropriate for the Cannonsville basin but must be verified empirically,
not assumed.

**Computational cost.** One SSI computation and one mean computation per trace
evaluation. SSI-3 requires a 3-month rolling sum and a per-month gamma fit
(applied once, not per evaluation). At the implemented evaluation costs for the SSI
pipeline, this constraint adds approximately 5-10% overhead per function evaluation.
This is acceptable.

**Constraint type.** Hard constraint for deviations beyond 0.30; soft constraint
between 0.15 and 0.30 (same two-tier structure as constraint 4.1).

### 4.4 Lag-1 autocorrelation tolerance reconciliation

**Current state.** The manuscript states tolerance 0.05; the code implements 0.30.
These must be reconciled before submission. The appropriate tolerance is neither value
as currently used.

**Calibration.** Bootstrap 1000 Kirsch traces of 30-year length (without optimization)
and compute the lag-1 autocorrelation of each monthly flow series. The 95% CI of this
distribution gives the sampling variability of lag-1 AC that is intrinsic to the
generator, not to the optimizer. Any constraint narrower than this CI would reject
statistically valid traces. Any constraint wider than this CI provides no additional
protection beyond what the generator already guarantees.

For monthly streamflow, the lag-1 autocorrelation of the full flattened monthly series
is dominated by the within-year seasonal cycle, which Kirsch preserves by construction.
The typical lag-1 AC of a 30-year monthly Kirsch trace is approximately 0.82-0.88 for
a basin with Cannonsville's seasonal amplitude. The 95% bootstrap CI of this statistic
at N=360 observations is approximately plus or minus 0.05-0.08. The 0.05 tolerance
stated in the manuscript is therefore approximately at the lower bound of what is
statistically defensible, and the 0.30 tolerance in the code is approximately six times
too wide. A calibrated tolerance in the range 0.08-0.12 (double to triple the bootstrap
standard error) is recommended as a defensible compromise between strict and permissive.

**Practical implication.** Because the within-year correlation structure is preserved
by Kirsch's Cholesky step, the lag-1 AC constraint is approximately redundant for any
trace that does not involve aggressive temporal clustering of low-flow year indices.
It becomes active only for pathological decision vectors. The constraint is worth
retaining as a catch-all but its calibration should be reported transparently in
Section 2.2.4.

### 4.5 Lag-2 autocorrelation constraint (optional)

**Rationale.** The lag-2 autocorrelation of the monthly series captures multi-month
persistence at the between-year transition that the lag-1 constraint and the Cholesky
step do not fully address. A trace with anomalous lag-2 AC could arise from the
optimizer selecting pairs of consecutive low-flow years without triggering the lag-1
constraint.

**Assessment.** The incremental protection provided by a lag-2 constraint beyond the
lag-1 constraint and the Hurst exponent constraint (Sections 4.2 and 4.4) is low.
The Hurst exponent integrates multi-scale autocorrelation structure more comprehensively
than any single lag pair. If the Hurst exponent constraint is implemented (Section 4.2),
the lag-2 constraint is unlikely to flag additional pathological traces and may be
omitted.

**Recommendation.** Skip if the Hurst exponent constraint is implemented. Add only
if the Hurst constraint is computationally prohibitive or statistically unstable at
30-year trace length.

---

## 5. Proposed Supporting Information Sensitivity Analysis

The purpose of this SI section is to provide empirical evidence that the constraint
set is neither overrestrictive (blocking large portions of the physically plausible
drought space) nor permissive (allowing statistically implausible traces into the
archive). The analysis structure follows the preemptive action requested in critique 2
of reviewer_critiques.md.

**Experiment structure.** Run four Cannonsville MOEA-FIND experiments, each at 100,000
function evaluations and five independent seeds:

- Baseline: current constraints (autocorr_tol=0.3, mean_tol=0.5, cv_tol=0.5, seasonal_tol=0.5)
- Proposed tight: annual mean tol=0.15, Hurst bootstrap CI, lag-1 AC bootstrap CI
- Proposed permissive: annual mean tol=0.30, Hurst 2x bootstrap CI, lag-1 AC 0.20
- No constraints: infeasibility disabled entirely

**Metrics reported for each experiment.**

(a) Archive size (number of non-dominated solutions retained). Smaller archive under
tighter constraints indicates the constraint is active and reducing the feasible
drought space.

(b) Nearest-neighbor coefficient of variation of the archive in drought characteristic
space. If the constraint is too tight, the archive may become sparser and less uniform.

(c) Fraction of function evaluations flagged as constraint-infeasible in each run.
This directly shows how often the constraints fire.

(d) Distribution of the constrained statistics (annual mean, Hurst H, lag-1 AC) across
archive members for each experiment. If the tight constraint is appropriate, this
distribution should be centered on the historical value with spread matching the
bootstrap CI.

(e) Median Hurst exponent of archive members vs. the historical Hurst exponent.
If the proposed constraint is correctly calibrated, these should agree within the
bootstrap CI.

**Table structure for SI.** One row per experiment, one column per metric. The table
explicitly shows whether the proposed constraints reduce archive quality (NN-CV increases)
or improve it (NN-CV unchanged but Hurst exponent and annual mean distributions are
closer to historical). The expected finding is that the proposed constraints have
negligible effect on the coverage metrics (because the optimizer naturally stays
within hydrologically plausible ranges) but eliminate the tail of implausible traces
that would otherwise appear at archive boundaries.

---

## 6. Priority Ranking

Constraints are ranked by the combination of (a) probability that a reviewer will
identify the gap as fatal, (b) ease of calibration from the existing data, and
(c) implementation complexity.

### Must-add before submission

**Priority 1: Fix the annual water balance constraint (Section 4.1).**
The current 50% tolerance on the annual mean is indefensible to any reviewer with
stochastic hydrology background. Replace with the bootstrap-calibrated 15% tolerance.
The code change is a single line in `compute_all_constraints`. The calibration
calculation requires one afternoon of bootstrap analysis on the historical record.
No reviewer in the synthetic hydrology community will accept a 50% excursion in annual
mean flow as "plausible." This is the highest-priority item.

**Priority 2: Implement the non-drought period constraint as stated in the manuscript
(Section 4.3).**
The manuscript claims a 15% tolerance on non-drought period flows, but the code does
not implement this. The discrepancy between the manuscript and the code is a submission-
blocking error independent of any methodological debate about the appropriate tolerance
value. This requires approximately 20 lines of code to implement correctly with the
event-detected drought/non-drought split.

**Priority 3: Reconcile the lag-1 autocorrelation tolerance (Section 4.4).**
The manuscript states 0.05 and the code uses 0.30. One of these values must be changed
to match the other, and the chosen value must be reported in the manuscript with a
sentence explaining how it was set. The bootstrap calibration described in Section 4.4
takes an afternoon and produces a defensible single-sentence justification for the
tolerance value. Without reconciliation, any reviewer who reads both the methods section
and the constraint source code will flag this as a methodological inconsistency.

### Strongly recommended

**Priority 4: Add a Hurst exponent constraint (Section 4.2).**
Wheeler et al. (2024) include the Hurst coefficient explicitly as an optimization
objective in a directly comparable framework. MOEA-FIND's failure to address inter-
annual persistence will be identified by any reviewer familiar with that paper or with
the stochastic hydrology literature on long-memory processes. The constraint is
computationally cheap (DFA on 360 observations) and the bootstrap calibration is
straightforward from the historical record and the Kirsch library. This is the most
theoretically significant gap in the current constraint set.

### Nice-to-have

**Priority 5: Annual flow quantile constraint (P10, P90).**
This provides additional protection against traces whose tail structure is distorted
even if the mean and CV are within bounds. It is useful but unlikely to be raised as
a fatal critique by itself, given that the annual mean and CV constraints already
bound the gross properties of the distribution. Add to the SI sensitivity analysis
but do not include in the main text unless the quantile constraint turns out to be
binding in experiments.

### Skip

**Lag-2 autocorrelation (Section 4.5).** Redundant given the Hurst exponent
constraint. Adds complexity without additional protection.

**Within-year monthly correlation structure.** Fixed by Kirsch's Cholesky step exactly.
Adding a constraint would be computationally redundant and would signal to reviewers
that the authors do not fully understand what their generator preserves.

**Monthly means and standard deviations.** Fixed by generator construction. Same
argument.

**Skewness.** Implicitly constrained by the log-normal framework and approximately
preserved by bootstrap resampling. Low probability of generating implausible skewness
values given the existing constraints on the mean and CV. Skip unless preliminary
experiments show skewness violations in the Cannonsville archive.

---

## Summary of Implementation Changes Required

The following changes are required in `src/constraints.py` and `manuscript_main_draft.md`
before submission:

1. In `compute_all_constraints`, replace `mean_tol=0.5` with a bootstrap-calibrated
   value (recommended 0.15) applied to the annual mean flow.
2. Add a `non_drought_period_constraint` function that (a) detects the drought period
   using SSI-3 event extraction applied to the synthetic trace, (b) computes the mean
   flow across non-drought months, and (c) returns a violation proportional to the
   fractional deviation from the historical non-drought mean. Apply this in
   `compute_all_constraints`.
3. Replace the lag-1 autocorrelation tolerance in `constraints.py` with the bootstrap-
   calibrated value and update the stated value in Section 2.2.4 of the manuscript.
4. Add a `hurst_exponent_constraint` function using DFA on the full 360-month series.
   Include this as a soft constraint (penalty) in `compute_all_constraints`.
5. Update Section 2.2.4 to accurately describe all four constraints and report their
   calibrated tolerances. The current text conflates the implemented constraints with
   the manuscript claims.
6. Add the constraint sensitivity analysis to Supporting Information as described in
   Section 5 of this document.

---

*References cited: Borgomeo et al. (2015, WRR); Brechmann et al. (2017, Stoch.
Environ. Res. Risk Assess.); Deb (2000, IEEE Trans. Evol. Comput.); Hurst (1951,
Trans. Amer. Soc. Civil Eng.); Kirsch et al. (2013, JWRPM); Montanari et al. (1997,
Stoch. Environ. Res. Risk Assess.); Nowak et al. (2010, WRR); Peng et al. (1994,
Phys. Rev. E); Salas (1993, in Maidment Handbook of Hydrology); Stedinger (1982,
WRR); Thomas and Fiering (1962); Wheeler et al. (2024, J. Hydrol. Eng.); Zaniolo
et al. (2023, Environ. Model. Softw.).*
