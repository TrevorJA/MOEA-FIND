# MOEA-FIND Introduction — Revised Draft (Round 2, post-critique-1)

*Branch: claude/zen-swanson. Round 1 critique applied 2026-04-14 (autonomous session).*
*Backup of pre-revision draft: introduction_revised_round1.md.*
*Round 1 critique: intro_critique_round1.md.*
*Last updated: 2026-04-14.*

---

## 1. Introduction

Multi-year meteorological drought is the binding supply-side risk for reservoir-
based water supply systems in water-stressed regions, and robust evaluation of
operating policies and drought management rules under deep uncertainty requires
synthetic streamflow ensembles that span a meaningful range of physically plausible
drought conditions. The many-objective robust decision making framework of Kasprzyk
et al. (2013), building on the robust decision making tradition of Lempert et al.
(2003) and the many-objective vulnerability analysis of Herman et al. (2015),
evaluates candidate operating policies by simulating each across an ensemble of
plausible future states of the world and aggregating system performance into
robustness measures that guide policy selection. The conclusions a planner can draw
about policy vulnerability depend directly on how thoroughly that ensemble spans the
drought hazard space defined by the severity, duration, seasonal timing, and
system-specific threshold characteristics that govern the operating rules and
robustness metrics for the basin under study. Bonham et al. (2025) demonstrate that
ensembles designed for coverage in the uncertain input space systematically
undersample challenging conditions in the downstream performance space, so that
input-space design criteria alone do not guarantee adequate coverage of the drought
hazard conditions that matter most for policy evaluation. The exploratory modelling,
synthetic hydrology, and scenario discovery communities have developed a portfolio of
complementary approaches to synthetic scenario design, each concentrating the
structured design effort at a different stage of the analysis pipeline, as shown in
Figure 1. This paper introduces a fourth approach, MOEA-FIND, that targets the
drought hazard outcome space directly through multi-objective evolutionary search
rather than as an emergent projection of a design applied upstream of that space.

The most widely applied approach to synthetic scenario design for water supply
vulnerability analysis is forward exploratory modeling (Figure 1, panel A), in
which the analyst applies a structured design in the uncertain input or parameter
space, propagates that design forward through a streamflow generator and simulation
model, and analyzes the resulting ensemble through downstream scenario discovery to
identify which combinations of input factors drive system failures. Moallemi et al.
(2020) formalize this workflow by distinguishing a decision space of candidate
operating policies, an uncertainty space of states of the world, and an outcome
space of system performance measures, where the structured sampling effort is
concentrated in the uncertainty space and the outcome space is populated by forward
propagation through the model chain. Latin hypercube sampling applied to uncertain
hydroclimate parameters, hidden Markov model hyperparameters, and human-system
demand multipliers generates ensembles whose members span the analyst-specified
uncertainty space (Hadjimichael et al., 2020; Quinn et al., 2018; Quinn et al.,
2019). The many-objective robust decision making formulation of Kasprzyk et al.
(2013) extends this framework to joint policy search and vulnerability assessment,
generating populations of non-dominated operating policies evaluated across the
scenario ensemble and applying robustness metrics that quantify performance
reliability across the uncertainty space. Downstream scenario discovery methods
including patient rule induction (Bryant and Lempert, 2010), classification and
regression tree analysis (Hadjimichael et al., 2020), and gradient-boosted tree
classifiers trained on the scenario ensemble identify the input factor combinations
that partition the ensemble into acceptable and unacceptable performance regions
(Quinn et al., 2018). Quinn et al. (2020) demonstrate that the scenario-neutrality
of input-space designs is undermined by the choice of input sampling range, which
conditions which states of the world the simulation encounters and which performance
outcomes are observed. For water supply systems under drought risk, the specific
limitation of the forward exploratory approach is that coverage of the drought
hazard space is an emergent projection of the input-space design through the
nonlinear generator map and drought event extraction procedure rather than a
directly controlled property of the ensemble, and no refinement of the input-space
design alone resolves this without coupling the design criterion to the drought
hazard space.

A second family of methods, broadly termed bottom-up vulnerability analysis or
scenario-neutral stress testing (Figure 1, panel B), addresses the design problem
from the question of where in the climate forcing space the system performance
crosses a management threshold rather than which input factors are uncertain. Brown
et al. (2012) developed the decision scaling framework, in which a systematic grid
is applied across the climate perturbation space, a reservoir simulation model maps
each perturbed forcing scenario to a system performance outcome, and the resulting
climate response function is analyzed to identify the threshold boundary in the
forcing space that separates acceptable from unacceptable performance. This design
makes the robustness assessment independent of any single climate projection, with
climate projections overlaid on the response function after the vulnerability
analysis is complete. Culley et al. (2016) extended the approach to the scenario-
neutral framework, in which a broad range of meteorological perturbations is applied
without prior conditioning on a climate model, generating a performance map that
spans the exposure space as broadly as the forcing design allows. Herman et al.
(2016) applied the scenario-neutral approach specifically to drought by developing
a weighted bootstrap resampling procedure that adjusts the block-selection
probabilities in a monthly nonparametric generator to produce synthetic traces with
increased drought frequency and severity relative to the unweighted historical
resampler. Their formulation characterizes drought events on the Standardized
Streamflow Index, evaluates the weighted ensemble against a reservoir simulation
model, and maps system performance across the space of drought frequency and
severity perturbations, with the generated traces remaining physically consistent
with the observed monthly block structure because the weight modification acts on
block-selection probabilities rather than on the values within each block. The
limitation of the bottom-up vulnerability approach from the perspective of the
present paper is that the structured design step is applied in the forcing variable
space, not in the drought hazard characteristic space. Drought severity, duration,
and timing appear in the generated ensemble as emergent consequences of the forcing
perturbation rather than as controlled design targets, and an analyst who wants to
understand which specific combinations of drought duration and severity are most
damaging must infer those combinations from the response function rather than
specifying them in advance as design inputs.

A third approach, library-and-subsample, generates a large pre-computed library of
synthetic traces from an unsteered stochastic generator and applies a structured
subsampling criterion to select a tractable ensemble for downstream simulation
(Figure 1, panel C). This approach shares the forward-propagation logic of panel A
but makes a qualitatively distinct methodological choice at the second stage,
applying a space-filling subsampling criterion within the coordinates of the pre-
generated library rather than designing the input space before generation. Bonham
et al. (2024) developed this approach for the Colorado River Basin by generating a
multi-scenario library of streamflow, demand, and initial storage combinations and
applying conditioned Latin hypercube sampling to select a subsample that covers the
generator input coordinates with near-uniform spacing. Their evaluation demonstrates
that input-space space-filling metrics are strong predictors of robustness ranking
preservation as subsample size decreases, so that a carefully selected small ensemble
recovers the ranking structure of the full library when the subsampling criterion is
well matched to the input space. This approach is computationally attractive because
the generator is run once to produce the library and the subsampling step is
inexpensive relative to regenerating traces or running the downstream simulation
model. The limitation that applies to library-and-subsample methods is that the
subsampling criterion is defined in the input coordinates of the generator rather
than in the drought hazard characteristic space the vulnerability analysis evaluates.
The generator map from input coordinates to drought characteristics is not in general
monotone or invertible, so a space-filling design in the input space is expected to
produce a non-uniform projection into the hazard space, oversampling the drought
conditions the generator produces most readily from the historical record and
undersampling the tails and combinations that are physically attainable but
statistically rare. Bonham et al. (2025) confirm this projection distortion as a
systematic bias in input-space ensembles across multiple case studies. The analyst
who wants structured coverage of the drought hazard space must either accept the
emergent projection of the input-space design or apply a secondary subsampling step
in the hazard coordinates, which requires a library large enough that the hazard
space is already well populated before subsampling.

A distinct and narrower question, not fully addressed by the three approaches above,
is whether synthetic scenario ensembles can be designed for structured coverage
directly in the drought hazard characteristic space, so that the analyst controls
which combinations of drought severity, duration, timing, and frequency appear in
the ensemble. Two prior lines of work address this question and together define the
gap that the present paper fills.

The first line applies bivariate and trivariate severity-duration-frequency analysis,
fitting a parametric joint distribution to the historical record of drought event
attributes and generating scenario ensembles by sampling from the fitted model.
Yevjevich (1967) introduced the run theory framework that defines drought duration,
intensity, and severity as measurable attributes of the interval during which a
streamflow index remains below a threshold, providing the event-level definitions
that the bivariate and trivariate severity-duration-frequency literature subsequently
adopted. Shiau (2006) showed that the joint distribution of drought duration and
severity extracted from a historical streamflow record can be represented by a
bivariate copula, and that copula-based sampling recovers the observed joint
exceedance probabilities of duration and severity more accurately than a bivariate
normal model with the same marginals. Serinaldi et al. (2009) extended the approach
to a trivariate copula connecting drought duration, mean deficit intensity, and peak
deficit intensity, and demonstrated that sampling from the fitted model produces
drought event sequences consistent with the historical joint frequency structure.
Salvadori and De Michele (2004) provided a unified treatment of bivariate copula
methods for joint extreme event frequency analysis that has been widely applied to
the joint distribution of drought event attributes in subsequent work. Copula-based
severity-duration-frequency methods are valid and well-validated for probabilistic
drought hazard characterization, and they provide a principled way to generate
scenarios with specified joint return periods of drought severity and duration. The
limitation for the design problem motivating the present paper applies specifically
to standalone severity-duration-frequency models: these methods sample from a
parametric model fitted to the observed event record and generate drought attribute
vectors directly, without coupling the scenario generation step to a bootstrap
streamflow generator. The generated scenarios describe joint exceedance probabilities
in the drought attribute space but do not produce physically realizable streamflow
time series whose monthly structure, seasonal cycle, and cross-month correlation
function are consistent with the historical record and with the generator used for
the non-drought portions of the ensemble. Hybrid approaches that embed copula
sampling within a conditional streamflow disaggregation model do exist and do not
share this limitation; the distinction motivating the present paper is therefore not
that all copula-based scenario generation is decoupled from the streamflow model,
but that standalone severity-duration-frequency models assume a parametric form for
the joint distribution of drought attributes and do not discover the achievable
region from the generator empirically.

The second line of work couples a streamflow generator directly to an optimization
procedure and produces one synthetic trace per user-specified drought target.
Borgomeo et al. (2015) introduced a simulated annealing procedure in which the
decision vector is a permutation of bootstrapped monthly flows and the objective
function combines deviations from a target vector of hydrological statistics with
penalties on non-drought-period flow statistics, producing one synthetic trace per
target that remains physically consistent with the observed record. The approach
demonstrated that directed search over the generator input space can produce
hydrologically realistic synthetic traces with prescribed statistical properties.
Zaniolo et al. (2024) adapted the procedure specifically to drought characteristic
targeting through the FIND method, defining three drought characteristics on the
Standardized Streamflow Index and minimizing a weighted objective function that
combines deviations from user-specified frequency, intensity, and duration targets
with secondary terms penalizing non-drought-period flow statistics, including mean
annual flow and flow percentiles. MOEA-FIND externalizes equivalent non-drought
plausibility controls into a separate constraint layer rather than objective terms,
which preserves the orthogonality of the drought characteristic objectives and
allows the Manhattan-distance construction to function; the design choice is made
explicit in Section 2. Wheeler et al. (in press) extended the single-objective
target-matching formulation to multi-site basins, expanding the target vector with
cross-site correlation matrices and producing one joint multi-site record per
high-dimensional target vector. The structural limitation common to all
target-matching methods is that each optimization run produces one synthetic trace
that matches one pre-specified target vector, so building an ensemble that spans a
range of drought characteristic combinations requires running the optimizer once per
target specification. The subset of the target combination space that is physically
attainable under the generator must be identified ex post by examining which
enumerated cells fail to converge, rather than being discovered automatically during
search. As the number of target drought characteristics grows, the enumeration grid
grows exponentially in the combination space while the fraction of the grid that is
physically attainable under the generator remains unknown until convergence is
attempted for each cell. The gap that motivates the present paper is therefore the
absence of a method that both discovers the feasible drought hazard region
automatically from an arbitrary stochastic generator and delivers structured coverage
of that region in a single optimization run, without requiring the analyst to
enumerate target combinations in advance and without assuming a parametric form for
the joint distribution of achievable drought characteristics.

This paper introduces MOEA-FIND, a fourth approach that places the structured design
step directly in the drought hazard characteristic space (Figure 1, panel D). The
method explores which drought signatures are achievable by recombination of observed
monthly values through the Kirsch-Nowak bootstrap generator (Kirsch et al., 2013;
Nowak et al., 2010) and does not extrapolate beyond the historical flow envelope.
Parametric extension for climate-change stress testing, through substitution of a
kappa-distribution generator or vine-copula model at the same coupling point, is
architecturally compatible and is deferred to follow-on work. Within that scope,
MOEA-FIND couples the Borg multi-objective evolutionary algorithm (Hadka and Reed,
2013) to the Kirsch-Nowak generator through a K+1-objective formulation in which a
single auxiliary Manhattan-distance objective causes every feasible synthetic trace
to be Pareto-optimal, not merely those on the non-dominated frontier of the drought
characteristic objectives. This formulation converts coverage of the feasible drought
hazard region from a post-hoc subsampling step into a property guaranteed by the
structure of the epsilon-dominance archive. The archive tiles the feasible region at
a resolution determined by the user-supplied epsilon vector, and convergence of the
archive constitutes the structured ensemble that is the deliverable of the method.
Two design properties require explicit statement. First, MOEA-FIND does not prescribe
what drought characteristics a basin can produce and does not redefine what a drought
is. The Kirsch-Nowak generator and its historical forcing remain the physical and
statistical model of drought emergence; the optimizer determines which combinations
of drought characteristics are achievable under that model without imposing a
parametric form on their joint distribution. Second, physical correlations among
drought characteristics, such as the tendency for longer events to accumulate larger
cumulative deficits and therefore to be more severe on average, are a feature of the
feasible drought hazard region rather than a limitation of the method. The method is
designed so that the Pareto-optimal archive adapts to the correlated structure of the
feasible region and covers its interior, including portions where physically coupled
characteristics co-occur; Section 3 evaluates empirically whether this coverage is
achieved on the Kirsch-Nowak generator applied to the Cannonsville basin. Marginal
and pairwise diagnostics for the archive members are reported in the coverage
analysis as a check on this property. The novelty of the present contribution is the
K+1-objective formulation itself: no existing method uses an auxiliary Manhattan-
distance objective that causes all feasible traces to be Pareto-optimal by
construction, and prior approaches to hazard-space targeting either fit parametric
distributions to the event record and sample from the fitted model, as in the
standalone copula-based severity-duration-frequency tradition, or optimize one
synthetic trace per user-specified target and discover the feasible region ex post,
as in the FIND tradition.

This paper makes four contributions to the literature on synthetic drought ensemble
design. The first is a reformulation of drought ensemble generation as a K+1-
objective evolutionary optimization problem in which an auxiliary Manhattan-distance
objective causes the epsilon-dominance archive of Borg MOEA to tile the full K-
dimensional feasible drought hazard region in one optimization run, automating the
feasibility discovery and coverage steps that require separate runs and ex-post
filtering in the predecessor literature. The second is the geometric argument that
justifies this behavior, showing that the auxiliary objective aligns all feasible
objective vectors onto a known codimension-one affine subset of the K+1-dimensional
objective space on which the entire feasible image is Pareto-optimal and the epsilon-
dominance archive tiles the subset at a user-specified resolution along each of the
K target drought characteristic axes. The third is an empirical verification of
interior-filling coverage on a constrained analytic benchmark across target
dimensionalities from two to six, ruling out shell-only and orthant-collapse failure
modes of the construction on a problem where the feasible region is convex and the
generator map is known; the analog verification on the Kirsch-Nowak application is
reported in Section 3. The fourth is a demonstration of the method on the
Cannonsville inflow with three custom drought characteristics, including a cyclic
peak timing feature, alongside a coverage comparison against equal-size library-
subsample baselines and a scenario discovery exercise conducted in drought hazard
space using operational failure labels derived from Pywr-DRB reservoir simulation as
the classification target. Section 2 defines the drought characterization pipeline,
the structured ensemble generation method, the analytic benchmark problem, the
Cannonsville case study, and the coverage and scenario discovery evaluation
protocols. Section 3 presents the three sets of empirical results. Section 4
discusses limitations, scope, and extensions. Section 5 concludes.

---

*End of introduction draft. Self-review notes follow below.*

---

## Self-review notes — round 1 critique responses (delete before submission)

**Critique 1 (gap statement mechanism):** Addressed. Para 7 now names the K+1
formulation explicitly: "no existing method uses an auxiliary Manhattan-distance
objective that causes all feasible traces to be Pareto-optimal by construction."
The novelty claim is now distinguishable from generic multi-objective drought
optimization.

**Critique 2 (SDF citation breadth):** Salvadori and De Michele (2004) added.
Hybrid SDF-generator approaches acknowledged with explicit scope statement: the
limitation applies to standalone severity-duration-frequency models, not to all
copula-based methods. A vine-copula drought citation was not added because no
specific paper could be verified with confidence; this remains an open item for
Trev.

**Critique 3 (FIND non-drought terms):** Addressed. Para 6 (FIND section) now
states that FIND's objective function includes both drought characteristic
deviations and secondary terms penalizing non-drought-period flow statistics, and
that MOEA-FIND externalizes equivalent controls into a constraint layer.

**Critique 4 (Wheeler et al. 2025 verifiability):** Changed to "Wheeler et al.
(in press)" to signal it is not yet in a database. Trev should verify the full
bibliographic entry and provide the DOI or preprint URL before submission.

**Critique 5 ("fourth approach" taxonomy):** Not restructured — Trev's Figure 1
spec explicitly requires four panels and the user specified this structure.
Instead, Para 4 now adds a sentence explaining why library-and-subsample is
methodologically distinct from panel A: "This approach shares the forward-
propagation logic of panel A but makes a qualitatively distinct methodological
choice at the second stage, applying a space-filling subsampling criterion within
the coordinates of the pre-generated library rather than designing the input space
before generation." This pre-empts the merge critique without violating the
four-panel structure Trev specified.

**Critique 6 (Bonham 2025 vs 2024 ambiguity):** Both Bonham et al. (2024) and
Bonham et al. (2025) PDF files exist in the literature/ directory, confirming
they are distinct publications. Para 4 uses (2024) for the Colorado River Basin
library-and-subsample work and (2025) for the projection-distortion finding.
Trev should confirm the full title and DOI of the 2025 paper before submission.

**Critique 7 (emergent property defusion reads as assertion):** Addressed. The
text now says "the method is designed so that the Pareto-optimal archive adapts
to the correlated structure... Section 3 evaluates empirically whether this
coverage is achieved." The distinction between design intent and demonstrated
result is maintained.

**Critique 8 (historical-envelope placement):** Addressed. The limitation now
appears as the second and third sentences of Para 7 (MOEA-FIND), immediately
after the introduction of the method name, before the technical mechanism is
described. It is no longer buried in a subordinate clause at the end.

**Critique 9 (opening motivation structure):** Partially addressed. Bonham et al.
(2025) is now cited in Para 1 as an empirical anchor for the claim that input-
space coverage is insufficient, rather than appearing only later in Para 2. The
opening sentence still starts with the broad motivation (drought risk) before
narrowing to the gap; a complete restructuring to lead with the empirical finding
would violate WRR stylistic expectations for an opening sentence.

**Critique 10 (Herman et al. 2016 Cholesky attribution):** Corrected. The
paragraph now says "block-selection probabilities rather than on the values
within each block," which accurately describes the Herman et al. weighted
bootstrap mechanism without attributing a Cholesky generator to a paper that
does not use one.

**Critique 11 (oversampling claim unsupported):** Softened. The text now says
"is expected to produce a non-uniform projection into the hazard space" rather
than asserting it as a demonstrated fact. The Bonham et al. (2025) citation
provides the empirical support for the specific bias claim.

**Critique 12 (Quinn et al. 2020 misattribution):** Corrected. Quinn et al.
(2020) now appears only in Para 2, where it is correctly attributed to the
claim that input-range conditioning undermines scenario-neutrality. It no longer
appears in Para 4 attributed to projection distortion from input space to
drought hazard space.

**Critique 13 (contributions list overstatement):** Addressed for the third
contribution, which now says "ruling out shell-only and orthant-collapse failure
modes of the construction on a problem where the feasible region is convex and
the generator map is known; the analog verification on the Kirsch-Nowak
application is reported in Section 3." This clearly limits what the analytic
benchmark proves. The fourth contribution's Pywr-DRB results are retained as a
contribution because the section3_3_redesign.md document specifies those
results; if they are not complete at submission, this item must be revised.

**Open items for Trev:**
1. Verify Wheeler et al. (in press) — provide full title, venue, DOI/preprint.
2. Verify Bonham et al. (2025) — provide full title and DOI (distinct from 2024).
3. Add one vine-copula drought frequency paper to the SDF paragraph if desired.
4. Confirm Pywr-DRB simulation results are available for the fourth contribution;
   if not, reframe as "workflow demonstration with proxy failure labels."
5. The "fourth approach" taxonomy (four panels) was retained per Trev's explicit
   specification. If reviewer insists on collapsing A and C, a one-paragraph
   response to review is already implicit in the Para 4 text.
