# Preemptive Peer-Reviewer Critique: MOEA-FIND

*Prepared 2026-04-14 against manuscript draft last updated 2026-04-14.
Purpose: anticipate the objections most likely from reviewers in the
synthetic hydrology, scenario discovery, and DMDU communities at
Water Resources Research and flag where the manuscript needs
strengthening before submission. Critiques are written in reviewer
voice. Four thematic sections cover methodological, conceptual,
statistical, and framing concerns. A closing "Hard Cases" section
identifies the five objections that most need authorial action.*

---

## Methodological Critiques

### 1. The formulation conflates two distinct purposes: it does not match a target and it does not guarantee coverage

**Critique.** The authors claim that the epsilon-dominance archive of
Borg MOEA "tiles" the feasible drought hazard region, but they never
prove that every epsilon-box that intersects the feasible image will
be populated given finite function evaluations. The analytic benchmark
demonstrates that the archive is interior-filling on a convex ball at
K=2 through K=6 under 30,000 function evaluations, but that is an
analytic benchmark with a two-to-six-dimensional decision space and a
known, convex, continuously differentiable feasible set. The real
Cannonsville case has 360 continuous decision variables, a generator
that is discontinuous in the decision-variable-to-drought-characteristic
map (the SSI event extraction involves threshold crossings), and a
feasible drought region of unknown shape and connectivity. There is no
convergence guarantee that scales from six dimensions to 360, and no
empirical evidence from the actual application is provided because the
Cannonsville HPC runs are pending. The manuscript cannot reasonably
claim structured coverage as a demonstrated result at this stage; it
can claim it only for the analytic benchmark.

**Reviewer would cite.** Reed et al. (2013, Advances in Water
Resources) demonstrate that MOEA performance degrades predictably with
decision-variable dimensionality and that many-objective convergence
diagnostics must be run on the specific problem, not inferred from
analytic benchmarks. Hadka and Reed (2015, Environmental Modelling and
Software) show that Borg requires problem-specific NFE calibration and
that restarts driven by epsilon-progress depend on the landscape having
accessible epsilon-boxes.

**Defensible response.** The authors are transparent about the
Cannonsville results being pending (HPC Phase beta), and the analytic
dimension sweep is a genuine contribution that rules out the failure
modes on the benchmark. However, the manuscript's abstract and key
points claim "demonstrated on the Cannonsville inflow" and present the
coverage comparison as Section 3.2, which is blank pending HPC. A
reviewer reading those sections will encounter pending placeholders for
both Sections 3.2 and 3.3, the two sections that contain the primary
empirical claim. This is a submission-readiness problem that is genuinely
hard to defend: the paper as drafted cannot be submitted without those
results. The appropriate action is to not submit until the HPC runs
complete.

**Preemptive action.** Before submission, ensure Sections 3.2 and 3.3
are fully populated with Cannonsville results, coverage metric tables,
and the scenario discovery comparison. Add an NFE and epsilon
sensitivity analysis for the Cannonsville problem specifically
(distinct from the analytic sweep), showing that the archive populates
the drought hazard region consistently across multiple random seeds and
does not collapse to a subset of the feasible space at 100,000 NFE.
Report multi-seed hypervolume convergence plots in Supporting
Information.

---

### 2. The soft-constraint implementation for plausibility is under-specified and potentially circular

**Critique.** Section 2.2.4 imposes two plausibility constraints:
lag-1 autocorrelation within 0.05 of the historical value and mean
non-drought annual flow within 15 percent of the historical mean. Both
are soft constraints implemented under the constraint-domination
criterion of Deb (2000). However, neither the tolerance values (0.05,
15%) nor the constraint penalty form are derived or calibrated. The
lag-1 autocorrelation constraint in particular may prevent the
optimizer from generating drought traces whose autocorrelation
structure differs from the historical record. Droughts with
anomalously persistent low-flow conditions (high positive
autocorrelation) or rapid recovery (low autocorrelation) are exactly
the cases that stress-testing practitioners need, and yet the
constraint as written prohibits them. Furthermore, the 15% mean-flow
tolerance is applied to the non-drought period, but the split between
drought and non-drought months changes as the optimizer generates more
extreme drought traces. It is not clear whether this constraint is
evaluated on a fixed set of months or on the event-detected split of
each trace.

**Reviewer would cite.** Borgomeo et al. (2015) include a
non-drought-period statistics penalty as an objective term in the
weighted-sum formulation, not as a hard constraint, explicitly to allow
the algorithm to trade off drought control against realism. Zaniolo et
al. (2024) similarly include autocorrelation and non-drought percentile
as objective terms with user-specifiable weights, and report that
weight retuning is necessary when targets move away from historical
conditions, implying that fixing these as constraints is even more
restrictive than the predecessor methods.

**Defensible response.** The constraints verify statistical plausibility
rather than identifying drought events: they check whether each synthetic
trace lies within the envelope of traces the Kirsch-Nowak generator could
produce naturally given an infinite ensemble, not whether it contains a
drought of a particular character. The lag-1 autocorrelation concern is
bounded by the Kirsch bootstrap, which enforces the historical correlation
structure by construction through the Cholesky factorisation: the
constraint is therefore approximately redundant for normal traces and
becomes active only for pathological decision vectors. The 15% non-drought
mean-flow tolerance is less well justified and should be calibrated against
the historical bootstrap distribution of non-drought flows. Tolerance
calibration is now run per experiment via
`workflows/diagnostics/diag_constraint_calibration.py` and reported in
the Supporting Information for each run; see DD-05 and DD-14 in
`governance/design_decisions.md`.

**Preemptive action.** Add a one-paragraph justification for both
tolerance values in Section 2.2.4 explaining how they were calibrated.
Show in Supporting Information that relaxing or tightening each
constraint by a factor of two does not materially change the coverage
metrics on the Cannonsville archive. Clarify whether the mean
non-drought flow constraint is computed on event-detected or fixed
calendar months.

---

### 3. The Manhattan-distance construction works only when the anti-ideal is placed correctly, and placing it correctly requires information not available before the run

**Critique.** The non-dominance argument in Section 2.2.3 requires
that the anti-ideal point D* satisfies D_j(x) <= D*_j for every
feasible x and every j. The manuscript states that D* is placed at the
historical Cannonsville drought characteristic vector. However, the
optimizer can and will generate synthetic traces whose drought
characteristics exceed the historical values in some dimensions (that
is the entire point of the method). If any archive member has
D_j(x) > D*_j, the absolute value in f_{K+1} does not reduce to
D*_j - D_j, the sum-to-constant identity fails, and the non-dominance
argument breaks down. The authors need to either place D* outside the
reachable image of the generator entirely, or demonstrate that the
historical vector is a provable upper bound on the feasible drought
characteristics. Neither is trivially satisfied: the bootstrap can
recombine months to produce droughts more severe than any observed
event.

**Reviewer would cite.** This is an internal consistency concern
rather than a literature conflict, but the affine identity of equation
(2) is the load-bearing claim that distinguishes MOEA-FIND from generic
multi-objective optimization of drought metrics. Any violation of the
anti-ideal placement assumption causes the method to degrade silently
into ordinary multi-objective search without the coverage guarantee.

**Defensible response.** The historical vector is used as D* in the
Cannonsville demonstration, which is a reasonable default, but the
manuscript should acknowledge that if the optimizer discovers traces
more extreme than historical in any characteristic, the construction
requires that D* be updated or placed conservatively above the
historical maximum with a buffer. The analytic benchmark uses D* = (3,
3, ..., 3) which is outside the K-ball by construction; an analogous
conservative placement for the hydrologic case is to set each D*_j to
the 99th percentile of the empirical bootstrap distribution over a
pilot library rather than the historical value. The residual diagnostic
(Section 3.1, standard deviation ~10^-16) verifies the identity on
the unconstrained analytic case, but an equivalent diagnostic should
be reported for the Cannonsville run.

**Preemptive action.** Add a paragraph in Section 2.2.3 explaining
the anti-ideal placement protocol for the hydrologic case. Report the
fraction of archive members that violate the anti-ideal placement
assumption D_j(x) <= D*_j in the Cannonsville run and show that it is
negligible. If any violations occur, explain how they are handled
(e.g., re-placement of D* at the maximum observed characteristic plus
a fixed buffer, with a re-run to verify the identity).

---

### 4. Convergence is assessed only by archive membership, not by hypervolume or gap-to-optimum

**Critique.** The manuscript validates coverage with nearest-neighbour
coefficient of variation, interior mass fraction, and orthant occupancy,
all of which are post-hoc diagnostics on the final archive. There is
no convergence diagnostic showing that the archive has stabilised: no
hypervolume-versus-NFE curves, no multi-seed variance in archive size
or coverage metrics, no comparison of archive quality between 50,000
and 100,000 NFE. Without these, I cannot assess whether 100,000 NFE is
sufficient for the 360-variable Cannonsville problem, whether the
archive is still growing at termination, or whether different random
seeds produce materially different archives.

**Reviewer would cite.** Reed et al. (2013) establish multi-seed
hypervolume convergence diagnostics as the standard for MOEA result
reporting in the water resources literature. Hadka and Reed (2015)
demonstrate that the number of function evaluations required to
achieve stationary hypervolume scales with decision-space dimension
and problem nonlinearity, and that 100,000 NFE is often insufficient
for high-dimensional environmental problems.

**Defensible response.** This is a standard methodological requirement
in MOEA papers and is completely defensible if addressed. The analytic
benchmark uses 30,000 NFE for K=6 (six-dimensional decision space),
and the paper notes that the Cannonsville run uses 100,000 NFE for
360 decision variables. The ratio of NFE to decision variables is
dramatically smaller for the hydrologic case, which is a genuine risk.

**Preemptive action.** Add Supporting Information figures showing
hypervolume-versus-NFE for the Cannonsville run across at least five
independent seeds. Report the mean and standard deviation of archive
size and coverage metrics across seeds. If the hypervolume has not
stabilised at 100,000 NFE, increase the budget or explicitly
acknowledge the limitation and frame the Cannonsville result as a
demonstration rather than a converged optimum.

---

### 5. Reproducibility is promised but the Borg dependency is a licensing barrier

**Critique.** The Open Research section states that the MOEA-FIND
implementation will be released on GitHub. However, Borg MOEA is a
licensed algorithm whose source code is not publicly distributable.
A researcher who wants to reproduce or extend the results will need to
independently obtain a Borg license, a non-trivial barrier. The
manuscript does not acknowledge this constraint, which contradicts the
implicit reproducibility claim.

**Reviewer would cite.** Hadka and Reed (2013) distribute Borg under
a research license that requires registration; the license is not GPL
or MIT. The WRR open data policy requires that reviewers can in
principle reproduce submitted results, which is not achievable if the
core optimizer is behind a license request.

**Defensible response.** This is a genuine limitation. The method is
algorithmically defined around epsilon-dominance archiving, and an
open alternative exists: the Python platypus library or the
MOEAFramework include epsilon-NSGA-II and epsilon-MOEA implementations
that use the same archiving mechanism and could substitute for Borg at
the cost of the adaptive operator selection. The manuscript should
acknowledge the Borg dependency, point users to the registration
process, and note that the wrapper is compatible with alternative
epsilon-dominance-based optimizers.

**Preemptive action.** Add a sentence in Section 2.2.2 noting that
Borg is available under a research license from the Reed group and
that the wrapper layer is designed to be optimizer-agnostic. Confirm
in the Open Research section that the repository provides an
example using a public epsilon-dominance optimizer so that reviewers
can run a toy version without a Borg license.

---

## Conceptual Critiques

### 6. Drought is an emergent system-level property, not a feature of streamflow alone; optimizing over streamflow drought metrics may not produce the drought scenarios that matter for a reservoir system

**Critique.** The method defines drought characteristics on the
Standardized Streamflow Index applied to the Cannonsville inflow. But
the operational drought that matters for reservoir management is
defined by reservoir storage, demand, and operating rules, not by
inflow alone. A sequence of moderate inflow deficits paired with
high demand or low initial storage may produce a catastrophic
operational drought, while a sequence of severe inflow deficits during
low-demand periods may cause no operational consequence. The method
covers the inflow drought hazard space uniformly, but it has no
mechanism to ensure that the covered region of inflow space maps to
a covered region of operational consequence space. Coverage of one
space is not coverage of the other when the mapping is nonlinear and
demand- and storage-dependent.

**Reviewer would cite.** Herman et al. (2016) make exactly this point:
"the framework for generating synthetic streamflow scenarios must be
sufficiently flexible to create extreme drought conditions not observed
in the historical record... to remain consistent with the bottom-up
approach, it must provide relatively simple parameters to adjust the
severity of scenarios in an exploratory manner." The emphasis is on
operational consequence, not on covering an inflow statistic space.
Hadjimichael et al. (2020, Earth's Future) frame scenario discovery
explicitly in terms of system performance outcomes for each
stakeholder, not in terms of input characteristics; their PRIM analysis
identifies factors conditional on performance failure, not on input
space uniformity. Bonham et al. (2025) note that "multiple studies
have shown this bias towards moderate SOWs can have complex,
difficult-to-anticipate impacts on subsequent factor mapping steps,"
specifically because input-space coverage does not guarantee
outcome-space coverage.

**Defensible response.** The manuscript is upfront that MOEA-FIND
covers the inflow drought hazard space by construction and that the
relationship between inflow drought characteristics and operational
consequences is mediated by the reservoir model and operating rules.
This is the same epistemic limitation shared by all streamflow-based
drought ensemble methods. The method's advantage over library
subsampling is that it covers the inflow drought hazard space
uniformly rather than density-weighted, which means the downstream
simulation samples the full range of inflow conditions rather than
oversampling mild conditions. But this is a partial answer: the
method cannot guarantee that the uniform inflow coverage translates
to informative operational coverage without the downstream Pywr-DRB
simulation results from the companion paper. This is a genuine
conceptual limitation that should be acknowledged honestly rather
than deflected.

**Preemptive action.** Add a paragraph in Section 4 explicitly
distinguishing inflow drought hazard space from operational consequence
space. State that MOEA-FIND provides structured coverage of the former
by construction and that whether this translates to improved
operational vulnerability discovery is an empirical question addressed
by the companion Pywr-DRB simulation paper. Frame this as a testable
hypothesis rather than a guaranteed outcome.

---

### 7. The method addresses coverage in drought hazard space, but it does not address what drought hazard space to cover; the anti-ideal choice is as consequential as the sampling method

**Critique.** The anti-ideal point D* determines the directions and
extent of the feasible drought hazard region that the optimizer
explores. In the Cannonsville case D* is placed at the historical
drought characteristic vector. This choice ensures exploration of
drought conditions at or below historical norms, which is appropriate
for validating against the observed record but may be too conservative
for stress testing. A reviewer doing bottom-up vulnerability analysis
wants to explore conditions beyond historical norms, not just
rearrangements of historical drought patterns. The method, as
described, is bounded by the historical record because the Kirsch
bootstrap samples only observed monthly values. The combination of
a historically grounded anti-ideal and a historically bounded generator
means the method explores a region that is by construction no more
extreme than what the historical record allows. The coverage
contribution is methodological precision within a bounded region, not
extension of the envelope.

**Reviewer would cite.** Borgomeo et al. (2015) explicitly note the
option to use perturbed target statistics to explore conditions beyond
the historical record. Herman et al. (2016) develop exactly this
extension, adjusting the resampling weights to increase drought
frequency and severity beyond historical levels. Zaniolo et al.
(2024) allow fractional perturbations of historical F, I, D targets
to explore non-stationary conditions. Quinn et al. (2020) argue that
scenario-neutral analyses are not neutral, specifically because
the range chosen for scenario generation conditions the vulnerabilities
discovered.

**Defensible response.** The manuscript acknowledges this in Section
4: the Kirsch bootstrap cannot generate drought signatures beyond the
historical envelope, and the parametric extension with kappa-vine
generators is identified as the path to extrapolation. This is an
honest framing. However, the Introduction presents the method as
addressing the gap of "which combinations are physically achievable
from the historical record and the underlying generator," which is a
correct statement of what the method delivers but may undersell the
generator's limitation relative to what practitioners need. The
contribution is mapping the achievable region under the Kirsch
generator, not extending it.

**Preemptive action.** Make the historical-envelope limitation
explicit in the abstract as well as in Section 4. Add a sentence
in the Introduction clarifying that the method is designed to
explore the set of drought signatures achievable by recombining
historical monthly values, not to extrapolate beyond the historical
record. Distinguish between the coverage contribution (sampling the
achievable region uniformly) and the generation contribution (defining
what is achievable).

---

### 8. The scenario discovery demonstration is circular: the failure label is defined on the same characteristics the optimizer covers

**Critique.** Section 3.3 defines a binary failure label as the joint
condition that mean drought severity exceeds the historical 80th
percentile and mean drought duration exceeds the historical 80th
percentile, then trains a gradient boosted tree classifier on this
label using the same (D_1, D_2) coordinates that MOEA-FIND is designed
to cover uniformly. It is tautological that a space-filling sample in
(D_1, D_2) produces a better-calibrated classifier of a threshold
defined in (D_1, D_2) than a density-weighted library subsample does.
The demonstration proves that space-filling sampling in the feature
space of a classifier yields better classification accuracy on a
threshold in that feature space, which is a property of the sampling
design, not a property of MOEA-FIND specifically. A rejection-based
LHS over the feasible region would produce the same result. The
demonstration does not show that MOEA-FIND's structured coverage
provides downstream value beyond what any feasible-region-restricted
space-filling sample would provide.

**Reviewer would cite.** Hadjimichael et al. (2020) and Quinn et al.
(2020) perform scenario discovery in parameter space because the
failure label is defined on system performance outcomes, not on the
parameter values used to generate the scenarios. The methodological
point is that the classifier must generalise from the training set to
new states of the world, which requires variation in the features that
matter for classification. Bonham et al. (2025) note that "multiple
studies have shown this bias towards moderate SOWs can have complex,
difficult-to-anticipate impacts on subsequent factor mapping steps,"
but their remediation is better sampling of the input space, not
sampling of the outcome space.

**Defensible response.** The circularity is real but is partially
addressed by framing the demonstration as a proof of concept for
hazard-space scenario discovery rather than a comparison of MOEA-FIND
to all possible samplers. The stronger demonstration would define the
failure label on a downstream operational metric (e.g., the fraction
of months in which Cannonsville reservoir storage falls below a
critical threshold in a Pywr-DRB simulation), and then compare
MOEA-FIND and library-LHS as training sets for a classifier of that
operational label. Without that downstream comparison, the scenario
discovery demonstration adds little to the coverage comparison of
Section 3.2. This is a genuine scope gap.

**Preemptive action.** Either strengthen the scenario discovery
demonstration to use an operational failure label from the Pywr-DRB
companion simulation, or be explicit that the Section 3.3
demonstration is a workflow illustration rather than an independent
validation of MOEA-FIND's decision-relevance. If the latter, reduce
the space given to Section 3.3 and reframe it accordingly. Acknowledge
in Section 4 that demonstrating improved downstream vulnerability
analysis using operational performance labels is deferred to the
companion paper.

---

### 9. The relationship to storyline-based scenario design is not addressed and is a missed opportunity to situate the contribution

**Critique.** The "physically self-consistent storyline" tradition
(Shepherd et al., 2018; cited in Hadjimichael et al. 2024) produces
scenarios that are physically plausible, interpretable narratives
rather than statistical samples. MOEA-FIND produces an archive of
synthetic traces that are indexed by drought characteristic vectors
but have no physical narrative structure. A practitioner using the
archive as an MORDM ensemble faces the question: which of these
archived traces should I present to stakeholders and how do I
communicate why each one is a relevant scenario? The method does not
address this communication problem, and a reviewer from the
storyline/scenario discovery community will note the gap.

**Reviewer would cite.** Hadjimichael et al. (2024, Earth's Future)
define scenario storylines as "physically self-consistent unfolding[s]
of plausible future events" and show that scenario discovery workflows
produce scenario descriptions that are more useful for multi-actor
planning than raw ensemble statistics. Moallemi et al. (2020b) identify
scenario generation as a named methodological fork (Stage II fork 3.2.2,
Framing future scenarios) and emphasise path dependency — the sensitivity
of downstream inferences to choices made at this fork. Herman et al.
(2015) characterise the same stage as Axis II (States of the World) in
their four-axis robustness taxonomy. Both frameworks position scenario
generation as load-bearing, and a reviewer from either tradition will
expect the manuscript to address interpretability as well as coverage.

**Defensible response.** MOEA-FIND produces a structured archive
whose members are indexed by interpretable drought characteristics,
which provides more narrative anchoring than a random library
subsample. A practitioner can select the archive member closest to
any particular storyline (e.g., "a drought similar to 1963-1965
but lasting one month longer") by proximity in drought characteristic
space. This is a feature of the method that the manuscript does
not exploit. The gap is real but manageable.

**Preemptive action.** Add a paragraph in Section 4 discussing the
relationship between the MOEA-FIND archive and storyline-based
scenario selection. Show that archive members can be selected by
proximity to historically named droughts in the drought characteristic
space, and that the method's structured coverage ensures that any
target storyline in the feasible region has a nearby archive member
within a distance bounded by the epsilon vector.

---

## Statistical Critiques

### 10. Drought characteristics are correlated by construction; the method may tile a region whose axes are not independent

**Critique.** Mean drought severity and mean drought duration are not
independent across synthetic traces generated by the Kirsch bootstrap.
Within an SSI event record, longer droughts tend to accumulate larger
cumulative deficits and therefore tend to be more severe. This physical
coupling means the achievable region in (D_1, D_2) space is not an
axis-aligned box but a non-convex region skewed along the
duration-severity correlation diagonal. MOEA-FIND will produce an
archive that covers the achievable region, but the coverage metrics
(nearest-neighbour CV, interior mass fraction) are computed against a
uniform distribution on the feasible region's bounding box. If the
feasible region is highly non-convex or elongated, the coverage metrics
will make the archive appear better-performing than it is relative to
what a practitioner needs.

**Reviewer would cite.** Zaniolo et al. (2024) note that FIND's three
drought characteristics are not fully independent and demonstrate
that the achievable grid of F-I-D combinations is a strict subset of
the rectangular target grid, with certain combinations infeasible under
their generator. Borgomeo et al. (2015) make the same observation about
the tradeoffs among streamflow statistics in the objective function
of the simulated annealing procedure.

**Defensible response.** The manuscript explicitly acknowledges that
the achievable region is a non-convex subset of the bounding box and
that MOEA-FIND reveals the feasible region as a byproduct of the
search. The coverage metrics are computed within the estimated feasible
region using the convex hull of the union of all samplers, which
is a reasonable approximation. However, a convex hull of a
non-convex feasible region over-estimates the feasible region and will
inflate coverage scores for both samplers, making the comparison
less discriminating.

**Preemptive action.** Add a diagnostic plot showing the marginal
distributions of D_1, D_2, and D_3 across the MOEA-FIND archive and
the library, with the cross-correlations among characteristics
annotated. Discuss the extent to which the correlation structure of
the archive mirrors the physical coupling expected from the SSI event
extraction, and whether the coverage metrics are sensitive to the
convex-hull approximation of the feasible region.

---

### 11. The correlation structure of the generated traces is checked only at lag 1; higher-order and multi-site dependence are not validated

**Critique.** The plausibility constraint in Section 2.2.4 checks
only lag-1 autocorrelation. The Kirsch-Nowak generator is explicitly
designed to reproduce the full 12-by-12 monthly correlation matrix
(Kirsch et al., 2013), but there is no validation that the
MOEA-FIND-steered generator preserves this structure. The optimizer
may satisfy the lag-1 constraint while producing traces with unrealistic
higher-lag persistence or unrealistic seasonal correlation patterns.
For a water supply system whose operating rules depend on seasonal
carryover storage, multi-month persistence that differs substantially
from the historical record can produce physically implausible drought
scenarios.

**Reviewer would cite.** Kirsch et al. (2013) establish reproduction
of the full monthly correlation matrix as the defining property of
their generator and demonstrate that it outperforms fractional Gaussian
noise specifically in this regard. Wheeler et al. (2025) include
monthly autocorrelation at lags 1 through 11 and interannual
autocorrelation as explicit terms in the objective function precisely
because relaxing them allows the optimizer to produce hydrologically
implausible sequences. Borgomeo et al. (2015) include monthly standard
deviation and multiple autocorrelation lags in the simulated annealing
target vector for the same reason.

**Defensible response.** The Kirsch bootstrap imposes the historical
correlation matrix by construction through the Cholesky factorisation
at every evaluation. The residual injection mode used in the main text
feeds Gaussian residuals back through the inverse normal score
transform, which preserves the Cholesky structure, so the lag-1
constraint may in fact be redundant. However, the manuscript does not
demonstrate this, and the claim depends on the injection mode being
implemented correctly.

**Preemptive action.** Add a Supporting Information figure showing
the full ACF and seasonal correlation matrix for a representative
sample of archive members compared to the historical record, following
the convention of Borgomeo et al. (2015, Supporting Information S1).
Clarify in Section 2.2.4 that the Cholesky factorisation enforces the
full monthly correlation structure by construction and that the lag-1
constraint is a secondary check rather than the primary mechanism for
correlation preservation.

---

### 12. Extrapolation beyond the historical envelope is presented as a limitation but the implications for the feasible-region shape are not analysed

**Critique.** Section 4 correctly notes that the Kirsch bootstrap
cannot generate droughts more severe than any recombination of
historical monthly values allows. The manuscript presents this as a
known limitation, but does not quantify how much of the "plausible
drought hazard space" is outside the historical envelope. If the
historically bounded feasible region is small relative to the region
that a practitioner considers plausible under climate change, then the
method is covering a small and potentially uninformative corner of the
space that actually matters. The comparison to library subsampling is
valid only within the shared feasible region; neither method accesses
the tails beyond the historical record.

**Reviewer would cite.** Herman et al. (2016) propose weighted
bootstrap specifically to extend drought severity and frequency beyond
historical levels, and demonstrate that the historical record for the
Research Triangle contains only a few severe drought events from which
policy-relevant tails cannot be reliably estimated. Bonham et al.
(2025) note that ensembles derived from historical resampling
"systematically underrepresent challenging conditions relative to
moderate ones" precisely because the historical distribution is
right-skewed toward moderate conditions.

**Defensible response.** The limitation is real and is shared by all
bootstrap-based methods, including FIND and the library baseline.
The method is not positioned as an extrapolation tool; it is
positioned as a structured sampler of the achievable region. The
honest framing is: MOEA-FIND discovers the boundary of what is
achievable from the historical record under the Kirsch generator and
covers it uniformly; the parametric extension with kappa-vine marginals
is needed for extrapolation and is deferred to follow-up work.

**Preemptive action.** Add a quantitative comparison in Section 4
or Supporting Information showing the extent of the historically
bounded feasible region relative to the expert-elicited plausible
range for the Cannonsville basin. If the feasible region is small
relative to the planning envelope, state this explicitly and
recommend the parametric extension for applications that require
coverage of the full planning envelope.

---

### 13. The coverage comparison requires the library baseline to be of adequate size, but ten thousand traces may not be sufficient for accurate feasible-region estimation in three dimensions

**Critique.** The Cannonsville comparison in Section 3.2 uses a
10,000-trace Kirsch library as the baseline. The feasible drought hazard
region is estimated as the convex hull of the union of all samplers.
If the 10,000-trace library has sparse coverage of the tails of the
(D_1, D_2, D_3) distribution, the convex hull estimate of the feasible
region will be under-estimated, and the coverage metrics will be
computed on a region that does not capture the true boundary. A
10,000-trace library in three dimensions with a block-bootstrap
generator that has limited tail coverage may underestimate the
feasible region by a significant margin.

**Reviewer would cite.** Bonham et al. (2024) use a 500-scenario
library for the Colorado River Basin and note that ensemble size
requirements depend on the dimensionality of the input space and the
type of robustness metric. Their framework is designed to find the
minimum adequate ensemble size, implying that there is no universal
answer. For a K=3 drought characteristic space with a non-convex
feasible region, 10,000 traces may or may not be sufficient; this
needs to be demonstrated.

**Defensible response.** Ten thousand traces is a reasonable starting
size for a K=3 coverage comparison, and the convex hull approximation
is a known limitation that is acknowledged in the manuscript. The
sensitivity of coverage metrics to library size should be shown
explicitly to bound the uncertainty in the feasible-region estimate.

**Preemptive action.** Add a library-size sensitivity analysis in
Supporting Information showing how the estimated feasible-region
boundary and the coverage metrics change as a function of library
size from 1,000 to 10,000 traces. Report the point at which the
boundary stabilises and confirm that 10,000 traces is in the stable
regime.

---

## Framing Critiques

### 14. The novelty claim that no prior method covers drought hazard space is stated forcefully but rests on an absence of evidence rather than a systematic review

**Critique.** The Introduction concludes: "No published method searches
or subsamples directly in drought hazard space for the purpose of
structured coverage across the target drought axes, and this is the
gap that the present work addresses." This is a strong claim. The
authors support it with reference to Hadjimichael et al. (2020),
Bonham et al. (2024), and Quinn et al. (2020), showing that those
papers operate in input or parameter space. However, the literature on
operational drought frequency analysis, stochastic drought generation,
and probabilistic scenario design is large, and the manuscript does
not include a systematic search to support the "no prior method" claim.
A single counterexample from outside the MORDM literature
(e.g., from the SDF curve or stochastic hydrology communities) could
undermine the novelty framing.

**Reviewer would cite.** The SDF (severity-duration-frequency) curve
literature (e.g., Yevjevich, 1967; Zelenhasic and Salvai, 1987) uses
bivariate or multivariate copula methods to characterise the joint
distribution of drought severity and duration, and some recent papers
(e.g., Shiau, 2006; Serinaldi et al., 2009) sample from that joint
distribution to generate scenario ensembles. Whether this constitutes
"structured coverage" in the sense used by the authors is debatable,
but it is a literature that should be engaged.

**Defensible response.** The SDF literature generates scenarios by
sampling from a parametric joint distribution of drought characteristics,
not by coupling a streamflow generator to a multi-objective optimizer
to discover the feasible region. The distinction is that MOEA-FIND
discovers the achievable region empirically from the generator rather
than assuming a parametric form for the distribution of drought
characteristics. This is a genuinely different approach. However, the
manuscript should acknowledge the SDF tradition and explain why
parametric sampling from a fitted bivariate distribution is not
sufficient for the use case.

**Preemptive action.** Add a brief paragraph in Section 1 engaging
the SDF and multivariate drought frequency analysis literature.
Explain that SDF-based sampling assumes a parametric joint distribution
fitted to historical events and does not account for the non-linear
mapping from streamflow generator decision variables to drought
characteristics. Position MOEA-FIND as complementary to SDF methods
rather than making a claim of absolute novelty.

---

### 15. The comparison with FIND is framed as "MOEA-FIND eliminates per-target weight retuning" but the manuscript does not demonstrate that weight retuning is actually a practical barrier

**Critique.** The Introduction argues that FIND requires per-target
weight retuning when targets move away from historical conditions,
citing Zaniolo et al. (2024). This is presented as a significant
operational limitation that MOEA-FIND overcomes. However, Zaniolo et
al. do not present weight retuning as a failure mode; they present it
as part of the calibration workflow. If weight retuning can be automated
by a one-time grid search at the beginning of a study, it is a minor
inconvenience, not a fundamental barrier. The manuscript would be
stronger if it showed empirically that FIND fails for some target
combination that MOEA-FIND succeeds on, rather than treating the
weight retuning requirement as a self-evident limitation.

**Reviewer would cite.** Zaniolo et al. (2024) report Experiment 3,
which enumerates a 5x5 grid of I-D combinations, each with independent
optimization runs. The per-regime weight retuning is mentioned in the
methods but the results show successful generation across the grid.
It is not obvious from the paper that weight retuning is prohibitively
costly for a systematic ensemble design.

**Defensible response.** The weight retuning limitation of FIND is
real but is a practical concern rather than a theoretical failure. The
stronger argument is that MOEA-FIND requires only one optimization
run (vs. 25 for the 5x5 grid) and returns the attainable subset as
a byproduct rather than requiring ex-post infeasibility detection.
The efficiency argument is more defensible than the "weight retuning
is impossible" argument.

**Preemptive action.** Reframe the FIND comparison in the Introduction
to emphasise computational efficiency and automatic feasibility
discovery rather than weight retuning difficulty. Compute the
approximate function-evaluation cost of running FIND for the same
coverage of the (D_1, D_2, D_3) space that MOEA-FIND provides, and
show that MOEA-FIND is more efficient even accounting for its higher
per-run NFE.

---

### 16. The Borg MOEA is presented as uniquely suited to the formulation, but this claim requires empirical support

**Critique.** Section 2.2.3 argues that the epsilon-dominance archive
of Borg MOEA is the specific mechanism that produces structured
coverage on the hyperplane, and that NSGA-III and MOEA/D are not
appropriate alternatives. However, NSGA-III could be configured with
reference directions aligned to the hyperplane, and epsilon-NSGA-II
from the platypus library uses the same epsilon-dominance archive as
Borg. The claim that "Borg MOEA rather than those alternatives" is
required is asserted but not demonstrated empirically. A reviewer
familiar with the MOEA literature will ask: have you run NSGA-III or
epsilon-NSGA-II on the same problem and shown that Borg produces
better coverage?

**Reviewer would cite.** Reed et al. (2013) compare Borg against
NSGA-II, GDE3, IBEA, MOEA/D, and others on environmental problems,
showing that Borg's adaptive operator selection and epsilon-progress
restart are the distinguishing mechanisms, not the epsilon-archive
alone. The epsilon-archive mechanism is shared with epsilon-NSGA-II.
The paper does not show that epsilon-NSGA-II, which also uses the
same archive, would fail where Borg succeeds on the MOEA-FIND problem.

**Defensible response.** The theoretical argument for Borg is sound:
epsilon-dominance archiving is the mechanism that tiles the hyperplane,
and Borg's adaptive operator selection and restart mechanisms are
advantages for the high-dimensional, non-convex Kirsch problem. The
claim that the epsilon-archive is sufficient while NSGA-III's
reference-direction mechanism is not is well-grounded in Section 2.2.3.
However, the claim would be strengthened by a brief empirical
comparison.

**Preemptive action.** The analytic benchmark results were produced
with EpsNSGAII (platypus) as a locally runnable stand-in for Borg MOEA.
When production Borg MOEA runs are completed on HPC for the same analytic
benchmark, the comparison between EpsNSGAII and Borg archive coverage
quality will be available implicitly. Add a Supporting Information note
reporting that EpsNSGAII and Borg MOEA produce equivalent interior-filling
coverage on the analytic benchmark because both use the same epsilon-
dominance archive mechanism, then show Borg's convergence advantage on
the higher-dimensional Cannonsville problem via multi-seed hypervolume
curves. This turns the uniqueness claim from an assertion into a supported
finding grounded in the specific epsilon-dominance archiving property.

---

### 17. The paper occupies an unusual position: it presents a methodology paper with most of its empirical content pending, targeting a high-impact applied journal

**Critique.** Water Resources Research is an applied hydrology journal.
The complete manuscript as submitted will need to demonstrate the
method on the Cannonsville application and show that it provides
practical value for the kind of reservoir management study that WRR
readers care about. The current draft has Sections 3.2, 3.3, and 5
as explicit pending placeholders. Sections 3.1 and 3.2 are described
in terms of "expected findings." A paper submitted with these
sections as pending placeholders will be desk-rejected. The paper's
contribution rests disproportionately on the analytic benchmark (the
K-ball), which proves the construction works in a synthetic test case.
That is a methods contribution more appropriate for a computational
intelligence venue than for WRR, unless the Cannonsville results are
present and convincing.

**Reviewer would cite.** No specific citation; this is an editorial
observation. But note that Borgomeo et al. (2015) and Zaniolo et al.
(2024) both present complete empirical results (Thames at Kingston
and Pit River in California respectively) before submission. Wheeler
et al. (2025) demonstrate results at 18 sites in the Eastern Nile.

**Defensible response.** This is a submission-readiness concern, not
a methodological weakness. The draft is clearly a working document
with explicit pending markers. It should not be submitted to WRR until
Sections 3.2, 3.3, and 5 are complete.

**Preemptive action.** Complete the HPC phases beta and gamma before
submission. Do not submit the manuscript to WRR with pending sections.

---

### 18. The framing as complementary to FIND rather than replacing it may be too gentle given that the two methods address the same problem

**Critique.** The Introduction is carefully written to position
MOEA-FIND as extending rather than replacing FIND, noting that FIND
"produces one synthetic trace per target triple" and that MOEA-FIND
"returns the grid of attainable targets as the deliverable of one run."
But the practical implication is that a practitioner choosing between
the two methods will use MOEA-FIND and not FIND: MOEA-FIND produces
a full ensemble with structured coverage in one run, while FIND
requires multiple runs and ex-post feasibility checking. The
complementary framing is diplomatically appropriate but may dilute
the novelty argument. A more direct framing would be: FIND is the
single-objective baseline that MOEA-FIND supersedes in the
specific use case of structured ensemble construction.

**Reviewer would cite.** Zaniolo et al. (2024) position FIND as the
first method to "directly and independently control drought frequency,
intensity, and duration," and that claim is not challenged by MOEA-FIND,
which does not control specific drought property values but discovers
the achievable set. The two methods have different goals and the
complementary framing is accurate.

**Defensible response.** The complementary framing is correct and
defensible. FIND is the right tool when a practitioner wants a trace
that matches a specific (F, I, D) target. MOEA-FIND is the right tool
when a practitioner wants an ensemble that covers the achievable (F, I, D)
space. The two use cases are genuinely different. The manuscript should
be clearer about when each method is appropriate.

**Preemptive action.** Add a decision table or paragraph in Section 4
specifying the conditions under which FIND is preferable to MOEA-FIND
and vice versa. This strengthens the complementary positioning and
makes the contribution more precise.

---

### 19. The decision-variable parameterisation is not compared against Zaniolo et al.'s block-replace-from-distribution operator, which is the direct predecessor

**Critique.** MOEA-FIND uses residual injection into the Kirsch
bootstrap as the decision-variable parameterisation. Zaniolo et al.
(2024) use a different operator: they replace a contiguous block of
months in the synthetic series with values drawn from a historical
n-month cumulative distribution function and disaggregated via k-NN.
This block-replacement operator is specifically designed to produce
coherent multi-month drought events and is validated against historical
drought statistics. The MOEA-FIND residual injection parameterisation
generates individual monthly residuals independently, which may produce
drought events that are temporally incoherent (e.g., alternating
drought and non-drought months within an event). The manuscript does
not compare the two parameterisations or justify why residual injection
is preferable for event-level drought targeting.

**Reviewer would cite.** Zaniolo et al. (2024, Methods Section 2.2)
describe the block-replacement operator and justify it as producing
temporally coherent drought events. The k-NN disaggregation step
ensures that replaced months transition smoothly from and to the
surrounding flow context.

**Defensible response.** The Kirsch bootstrap with Cholesky
correlation enforcement naturally produces temporally coherent
sequences without requiring a block-replacement operator, because the
correlation structure is imposed across all months simultaneously.
The lag-1 autocorrelation constraint provides an additional check.
However, the manuscript should show empirically that the
MOEA-FIND-generated drought events are temporally coherent by
reporting the distribution of within-event SSI trajectories across
archive members.

**Preemptive action.** Add a Supporting Information figure showing
SSI-3 time series for a sample of archive members spanning the range
of (D_1, D_2, D_3) values, demonstrating that the events are
temporally coherent and that drought onset and recovery transitions
are physically plausible. Compare the temporal structure of the
generated events to historical Cannonsville drought events.

---

## Hard Cases

The following five objections are the most difficult for the authors
to defend without significant scope changes or acknowledgment of
limitations that are currently under-stated. They are listed in
decreasing order of severity.

### HC-1. Sections 3.2, 3.3, and 5 are pending; the paper cannot be submitted to WRR in its current form

This is the most serious concern and is not methodological but
editorial. The manuscript's two primary empirical contributions
(Cannonsville coverage comparison and scenario discovery demonstration)
are explicitly marked as pending HPC completion. A WRR submission
requires complete results. The paper should not be submitted until
those sections are filled with actual numbers. In the meantime, the
analytic benchmark is a legitimate partial contribution, but it is
insufficient for a WRR paper without the hydrologic application. No
response can address this without completing the HPC runs.

### HC-2. The scenario discovery demonstration is circular and adds no independent validation of the method's decision-relevance

> **RESOLVED 2026-04-14.** Section 3.3 has been redesigned to use a
> Pywr-DRB-derived operational failure label (Hashimoto reliability
> and vulnerability against an FFMP-aligned threshold), so the failure
> label is no longer defined on the same coordinate system that
> MOEA-FIND optimises over. Production protocol and references are in
> `planning/experiment_plan.md` Part A.

The original concern (retained for context): the failure label in Section
3.3 was originally defined on the same (D_1, D_2) coordinates that
MOEA-FIND covers uniformly. Showing that a space-filling sample in a
coordinate system is a better training set for a classifier of a threshold
in that coordinate system proves nothing about the method's value to
practitioners. The demonstration needs an operational failure label
(reservoir performance from a Pywr-DRB simulation) to be informative.

### HC-3. The emergent-drought critique (Critique 6) exposes a gap between inflow coverage and operational consequence coverage that is not bridged in this paper

The method covers the inflow drought hazard space uniformly. Whether
this translates into improved operational vulnerability discovery
depends on the downstream simulation model and operating rules, and
that question is deferred to the companion Pywr-DRB paper. The WRR
submission needs to either include the companion results or explicitly
acknowledge that the operational relevance claim is untested. A methods
paper that claims practical value for reservoir management but does not
demonstrate reservoir management outcomes is vulnerable to a
"so what" rejection from WRR editors.

### HC-4. The anti-ideal placement assumption (Critique 3) is not verified for the hydrologic case and could silently invalidate the non-dominance argument

If the Kirsch bootstrap generates any trace whose drought
characteristics exceed the historical anti-ideal in any dimension,
the sum-to-constant identity fails and the Pareto-optimality argument
breaks down. The paper verifies the identity for the analytic benchmark
to machine precision but does not report the analogous diagnostic for
Cannonsville. Because the bootstrap can recombine monthly values to
produce droughts more extreme than any historical event in some
characteristics, this is not merely a theoretical concern. The authors
must either report the diagnostic for Cannonsville or modify the
anti-ideal placement protocol to guarantee that no feasible trace
exceeds D* in any characteristic.

### HC-5. The computational cost comparison against FIND and library-subsampling is not quantified

The manuscript presents MOEA-FIND as more efficient than running FIND
25 times for a 5x5 grid of targets, but it does not report the actual
wall-clock or NFE cost of the Cannonsville MOEA-FIND run or compare
it to the cost of generating the 10,000-trace library and subsampling
it. If MOEA-FIND requires 500,000 NFE (the number is not specified for
the Cannonsville case) and each evaluation takes several seconds, the
total cost may be comparable to or exceed the cost of the library
baseline. Reviewers and practitioners need this comparison to assess
the practical advantage of the method. Without it, the efficiency
claim is unsubstantiated.
