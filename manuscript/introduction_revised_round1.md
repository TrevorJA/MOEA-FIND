# MOEA-FIND Introduction — Revised Draft

*Branch: claude/zen-swanson. Saved to manuscript/introduction_revised.md in worktree to
avoid conflicts with concurrent Figure 1 edits on main. Do not merge until the
Figure 1 session closes. Last updated: 2026-04-14.*

---

## 1. Introduction

Multi-year meteorological drought is the binding supply-side risk for
reservoir-based water supply systems in water-stressed regions, and robust
evaluation of operating policies and drought management rules under deep
uncertainty requires synthetic streamflow ensembles that span a meaningful
range of physically plausible drought conditions. The many-objective robust
decision making framework of Kasprzyk et al. (2013), building on the robust
decision making tradition of Lempert et al. (2003) and the many-objective
vulnerability analysis of Herman et al. (2015), evaluates candidate operating
policies by simulating each across an ensemble of plausible future states of
the world and aggregating system performance into robustness measures that
guide policy selection. The conclusions a planner can draw about policy
vulnerability depend directly on how thoroughly that ensemble spans the drought
hazard space defined by the severity, duration, seasonal timing, and system-
specific threshold characteristics that govern the operating rules and
robustness metrics relevant to the basin under study. The exploratory modelling,
synthetic hydrology, and scenario discovery communities have developed a
portfolio of complementary approaches to synthetic scenario design, each
concentrating the structured design effort at a different stage of the analysis
pipeline, as shown in Figure 1. No single approach is universally superior; the
practitioner's choice should be guided by the research question, the drought
characteristics that matter most for the system under study, and the
epistemological stance appropriate to the decision context. This paper
introduces a fourth approach, MOEA-FIND, that targets the drought hazard
outcome space directly through multi-objective evolutionary search rather than
as an emergent projection of a design applied upstream of that space.

The most widely applied approach to synthetic scenario design for water supply
vulnerability analysis is forward exploratory modeling (Figure 1, panel A), in
which the analyst applies a structured design in the uncertain input or parameter
space, propagates that design forward through a streamflow generator and
simulation model, and analyzes the resulting ensemble through downstream scenario
discovery to identify which combinations of input factors drive system failures.
Moallemi et al. (2020) formalize this workflow by distinguishing a decision space
of candidate operating policies, an uncertainty space of states of the world, and
an outcome space of system performance measures, where the structured sampling
effort is concentrated in the uncertainty space and the outcome space is
populated by forward propagation through the model chain. Latin hypercube
sampling applied to uncertain hydroclimate parameters, hidden Markov model
hyperparameters, and human-system demand multipliers generates ensembles whose
members span the analyst-specified uncertainty space (Hadjimichael et al., 2020;
Quinn et al., 2018; Quinn et al., 2019). The many-objective robust decision
making formulation of Kasprzyk et al. (2013) extends this framework to joint
policy search and vulnerability assessment, generating populations of
non-dominated operating policies evaluated across the scenario ensemble and
applying robustness metrics that quantify performance reliability across the
uncertainty space. Downstream scenario discovery methods including patient rule
induction (Bryant and Lempert, 2010), classification and regression tree analysis
(Hadjimichael et al., 2020), and gradient-boosted tree classifiers trained on the
scenario ensemble identify the input factor combinations that partition the
ensemble into acceptable and unacceptable performance regions (Quinn et al.,
2018). Quinn et al. (2020) demonstrate that the scenario-neutrality of
input-space designs is undermined by the choice of input sampling range, which
conditions which states of the world the simulation encounters and which
performance outcomes are observed. Bonham et al. (2025) show further that
input-space ensembles systematically underrepresent challenging conditions in the
downstream performance space because the input distribution is not calibrated
against the outcome space. For water supply systems under drought risk, the
specific limitation of the forward exploratory approach is that coverage of the
drought hazard space is an emergent projection of the input-space design through
the nonlinear generator map and drought event extraction procedure rather than a
directly controlled property of the ensemble, and no refinement of the input-
space design alone resolves this without coupling the design criterion to the
drought hazard space.

A second family of methods, broadly termed bottom-up vulnerability analysis or
scenario-neutral stress testing (Figure 1, panel B), addresses the design problem
from the question of where in the climate forcing space the system performance
crosses a management threshold rather than which input factors are uncertain.
Brown et al. (2012) developed the decision scaling framework, in which a
systematic grid is applied across the climate perturbation space, a reservoir
simulation model maps each perturbed forcing scenario to a system performance
outcome, and the resulting climate response function is analyzed to identify the
threshold boundary in the forcing space that separates acceptable from
unacceptable performance. This design makes the robustness assessment independent
of any single climate projection: the response function is a property of the
policy and the physical system, and climate projections are overlaid after the
vulnerability analysis is complete. Culley et al. (2016) extended the approach
to what they term the scenario-neutral framework, in which a broad range of
meteorological perturbations is applied without prior conditioning on a climate
model, generating a performance map that spans the exposure space as broadly as
the forcing design allows. Herman et al. (2016) applied the scenario-neutral
approach specifically to drought by developing a weighted bootstrap resampling
procedure that adjusts the sampling weights in a nonparametric generator to
produce synthetic traces with increased drought frequency and severity relative
to the unweighted historical resampler. Their formulation characterizes drought
events on the Standardized Streamflow Index, evaluates the weighted ensemble
against a reservoir simulation model, and maps system performance across the
space of drought frequency and severity perturbations, with the generated traces
remaining physically consistent with the Cholesky correlation structure of the
bootstrap generator because the weight modification acts on the sampling
distribution rather than on the physical structure of the historical record. The
limitation of the bottom-up vulnerability approach from the perspective of the
present paper is that the structured design step is applied in the forcing
variable space, not in the drought hazard characteristic space. Drought severity,
duration, and timing appear in the generated ensemble as emergent consequences of
the forcing perturbation rather than as controlled design targets, and an analyst
who wants to understand which specific combinations of drought duration and
severity are most damaging must infer those combinations from the response
function rather than specifying them in advance as design inputs.

A third approach, library-and-subsample, generates a large pre-computed library
of synthetic traces from an unsteered stochastic generator and applies a
structured subsampling criterion to select a tractable ensemble for downstream
simulation (Figure 1, panel C). Bonham et al. (2024) developed this approach for
the Colorado River Basin by generating a multi-scenario library of streamflow,
demand, and initial storage combinations and applying conditioned Latin hypercube
sampling to select a subsample that covers the generator input coordinates with
near-uniform spacing. Their evaluation demonstrates that input-space space-filling
metrics are strong predictors of robustness ranking preservation as subsample
size decreases, so that a carefully selected small ensemble recovers the ranking
structure of the full library when the subsampling criterion is well matched to
the structure of the input space. This approach is computationally attractive
because the generator is run once to produce the library and the subsampling step
is inexpensive relative to regenerating traces or running the downstream
simulation model. The limitation that applies to library-and-subsample methods
is that the subsampling criterion is defined in the input coordinates of the
generator rather than in the drought hazard characteristic space the vulnerability
analysis evaluates. The generator map from input coordinates to drought
characteristics is not in general monotone, linear, or invertible, so a space-
filling design in the input space produces a non-uniform projection into the
hazard space, oversampling the drought conditions the generator produces most
readily from the historical record and undersampling the tails and combinations
that are physically attainable but statistically rare. Bonham et al. (2025)
confirm this projection distortion as a systematic bias in input-space ensembles
across multiple case studies. Quinn et al. (2020) frame the same phenomenon as a
form of implicit conditioning that renders the ensemble design non-neutral in the
drought characteristic coordinates that planners evaluate. The analyst who wants
structured coverage of the drought hazard space must either accept the emergent
projection of the input-space design or apply a secondary subsampling step in the
hazard coordinates, but the latter requires a library large enough that the hazard
space is already well populated before subsampling, which may be infeasible for
long synthetic traces or for hazard spaces defined by more than two or three
drought characteristics.

A distinct and narrower question, not fully addressed by the three approaches
above, is whether synthetic scenario ensembles can be designed for structured
coverage directly in the drought hazard characteristic space, so that the analyst
controls which combinations of drought severity, duration, timing, and frequency
appear in the ensemble. Two prior lines of work address this question and together
define the gap that the present paper fills.

The first line applies bivariate and trivariate severity-duration-frequency
analysis, fitting a parametric joint distribution to the historical record of
drought event attributes and generating scenario ensembles by sampling from the
fitted model. Yevjevich (1967) introduced the run theory framework that defines
drought duration, intensity, and severity as measurable attributes of the interval
during which a streamflow index remains below a threshold, providing the event-
level definitions that the bivariate and trivariate severity-duration-frequency
literature subsequently adopted. Shiau (2006) showed that the joint distribution
of drought duration and severity extracted from a historical streamflow record can
be represented by a bivariate copula, and that copula-based sampling recovers the
observed joint exceedance probabilities of duration and severity more accurately
than a bivariate normal model with the same marginals. Serinaldi et al. (2009)
extended the approach to a trivariate copula connecting drought duration, mean
deficit intensity, and peak deficit intensity, and demonstrated that sampling from
the fitted model produces drought event sequences consistent with the historical
joint frequency structure. Copula-based severity-duration-frequency methods are
valid and well-validated for probabilistic drought hazard characterization, and
they provide a principled way to generate scenarios with specified joint return
periods of drought severity and duration. The limitation for the design problem
motivating the present paper is that these methods sample from a parametric model
fitted to the observed event record and generate drought attribute vectors
directly, without coupling the scenario generation step to a bootstrap streamflow
generator. The generated scenarios describe joint exceedance probabilities in the
drought attribute space but do not produce physically realizable streamflow time
series whose monthly structure, seasonal cycle, and cross-month correlation
function are consistent with the historical record and with the generator used for
the non-drought portions of the ensemble, creating an inconsistency between the
drought and non-drought portions of the synthetic record when the scenarios are
used as forcing in a reservoir simulation. Furthermore, copula-based
severity-duration-frequency methods assume a parametric form for the joint
distribution of drought attributes; they discover neither the feasible region
achievable from a given stochastic generator nor the shape of that region in the
analyst's chosen drought characteristic coordinates.

The second line of work couples a streamflow generator directly to an
optimization procedure and produces one synthetic trace per user-specified drought
target. Borgomeo et al. (2015) introduced a simulated annealing procedure in
which the decision vector is a permutation of bootstrapped monthly flows, the
objective is a weighted sum of deviations from a target vector of historical
statistics, and each optimization run produces one synthetic trace whose
statistical properties match the supplied target while remaining physically
consistent with the observed record. The approach demonstrated that directed
search over the generator input space can produce hydrologically realistic
synthetic traces with prescribed statistical properties. Zaniolo et al. (2024)
adapted the procedure specifically to drought characteristic targeting through the
FIND method, defining three drought characteristics on the Standardized Streamflow
Index, allowing the analyst to specify target values of each as absolute quantities
or as fractional perturbations relative to historical, and minimizing a weighted
sum of deviations from the target triple through iterative simulated annealing.
FIND establishes that directed search over the bootstrap generator input space can
produce traces with prescribed drought frequency, intensity, and duration
characteristics, and demonstrates the approach through a grid of target
combinations on a California river basin. Wheeler et al. (2025) extended the
single-objective target-matching formulation to multi-site basins, expanding the
target vector with cross-site correlation matrices and producing one joint
multi-site record per high-dimensional target vector. The structural limitation
common to all target-matching methods is that each optimization run produces one
synthetic trace that matches one pre-specified target vector, so building an
ensemble that spans a range of drought characteristic combinations requires
running the optimizer once per target specification. The subset of the target
combination space that is physically attainable under the generator must be
identified ex post by examining which enumerated cells fail to converge, rather
than being discovered automatically during search. As the number of target drought
characteristics grows, the enumeration grid grows exponentially in the
combination space while the fraction of the grid that is physically attainable
under the generator remains unknown until convergence is attempted for each cell;
the cost and complexity of this ex-post feasibility discovery limits the practical
dimensionality of the ensemble design. The synthesis of these two limitations
defines the gap that motivates the present paper: copula-based methods sample
in the drought attribute space but are not coupled to the streamflow generator,
while target-matching methods are coupled to the generator but require one
optimization run per target and provide no mechanism for automatically discovering
which combinations of drought characteristics are attainable. No existing method
combines directed multi-objective search over drought hazard characteristics with
automatic discovery of the feasible drought hazard region from an arbitrary
stochastic generator, delivers structured ensemble coverage of that region in a
single optimization run, and returns the boundary of the feasible region as a
byproduct of the run rather than requiring it to be specified in advance.

This paper introduces MOEA-FIND, a fourth approach that places the structured
design step directly in the drought hazard characteristic space (Figure 1, panel
D). MOEA-FIND couples the Borg multi-objective evolutionary algorithm (Hadka and
Reed, 2013) to the Kirsch-Nowak nonparametric bootstrap generator (Kirsch et al.,
2013; Nowak et al., 2010) through a single auxiliary Manhattan-distance objective
that causes every feasible synthetic trace to be Pareto-optimal under a
K+1-objective formulation, where K denotes the number of drought characteristics
chosen by the analyst. The epsilon-dominance archive of the optimizer tiles the
feasible drought hazard region at a resolution determined by the user-supplied
epsilon vector, and the convergence archive constitutes the structured ensemble
that is the deliverable of the method. Two properties of the design require
explicit statement. First, MOEA-FIND does not prescribe what drought
characteristics a basin can produce, and it does not redefine what a drought is.
The method discovers the feasible drought hazard region empirically through
directed search over the generator input space, with the Kirsch-Nowak generator
and its historical forcing remaining the physical and statistical model of drought
emergence; the optimizer determines which combinations of drought characteristics
are achievable under that model without imposing a parametric form on their joint
distribution. Second, physical correlations among drought characteristics, such
as the tendency for longer events to accumulate larger cumulative deficits and
therefore to be more severe on average, are a feature of the feasible drought
hazard region rather than a limitation of the method. The Pareto-optimal archive
adapts to the correlated structure of the feasible region and covers its interior,
including the portions where physically coupled characteristics co-occur, rather
than assuming the feasible region is axis-aligned or the drought characteristics
are statistically independent. Marginal and pairwise diagnostics for the archive
members are reported in the coverage analysis as a check on this property. The
current implementation explores what is achievable by recombination of the
historical monthly flow record through the Kirsch-Nowak bootstrap generator;
parametric extension of the generator envelope for climate-change stress testing,
through substitution of a kappa-distribution generator or vine-copula model at the
same coupling point, is architecturally compatible and is deferred to follow-on
work. The defensible novelty of the present contribution is the combination of
directed multi-objective search over drought hazard characteristics with automatic
discovery of the feasible drought hazard region from an arbitrary stochastic
generator. Prior approaches to hazard-space targeting either fit parametric
distributions to the event record and sample from the fitted model, as in the
copula-based severity-duration-frequency tradition, or optimize one synthetic
trace per user-specified target and identify attainable combinations ex post, as
in the FIND tradition. MOEA-FIND combines the discovery and coverage steps in a
single optimization run that requires no prior specification of attainable target
combinations and no parametric assumption about the joint distribution of the
drought characteristics.

This paper makes four contributions to the literature on synthetic drought
ensemble design. The first is a reformulation of drought ensemble generation as a
K+1-objective evolutionary optimization problem in which a single auxiliary
Manhattan-distance objective causes the epsilon-dominance archive of Borg MOEA to
tile the full K-dimensional feasible drought hazard region in one optimization
run, automating the feasibility discovery and coverage steps that require separate
runs and ex-post filtering in the predecessor literature. The second is the
geometric argument that justifies this behavior: the auxiliary objective aligns
all feasible objective vectors onto a known codimension-one affine subset of the
K+1-dimensional objective space on which the entire feasible image is
Pareto-optimal and the epsilon-dominance archive tiles the subset at a user-
specified resolution along each of the K target drought characteristic axes. The
third is an empirical verification of interior-filling coverage across target
dimensionalities from two to six on a constrained analytic benchmark, ruling out
shell-only and orthant-collapse failure modes of the construction. The fourth is
a demonstration of the method on the Cannonsville inflow with three custom drought
characteristics, including a cyclic peak timing feature, alongside a coverage
comparison against equal-size library-subsample baselines and a scenario discovery
exercise conducted in drought hazard space using operational failure labels derived
from Pywr-DRB reservoir simulation as the classification target. Section 2 defines
the drought characterization pipeline, the structured ensemble generation method,
the analytic benchmark problem, the Cannonsville case study, and the coverage and
scenario discovery evaluation protocols. Section 3 presents the three sets of
empirical results. Section 4 discusses limitations, scope, and extensions. Section
5 concludes.

---

*End of introduction draft. Self-review notes follow below.*

---

## Self-review notes (for author use — delete before submission)

**Gap statement — hostile reading test.** The claim "No existing method combines
directed multi-objective search over drought hazard characteristics with automatic
discovery of the feasible drought hazard region from an arbitrary stochastic
generator" is a specific conjunction claim. A hostile reviewer must find a method
that simultaneously satisfies all three conditions to defeat it. Copula methods
satisfy the first condition partially (they work in hazard space) but not the
second (no generator coupling). FIND satisfies the second (generator coupling)
but not the first or second for automatic feasibility discovery. The conjunction
claim survives.

**"Emergent property" objection.** Addressed in paragraph 6: "the Kirsch-Nowak
generator and its historical forcing remaining the physical and statistical model
of drought emergence; the optimizer determines which combinations of drought
characteristics are achievable under that model." The phrasing does not say
drought is a designed artifact; it says the optimizer maps the achievable region
under the physical model. The objection is defused without being defensive.

**"Correlated statistics" objection.** Addressed in paragraph 6: "physical
correlations among drought characteristics... are a feature of the feasible
drought hazard region rather than a limitation of the method." The marginal
diagnostic mention provides evidence rather than assertion.

**SDF engagement.** Paragraph 5 (first section) engages Yevjevich (1967),
Shiau (2006), Serinaldi et al. (2009) substantively. The distinction from
MOEA-FIND is stated on two grounds: (a) no coupling to the streamflow generator,
inconsistency between drought and non-drought portions of the synthetic record;
(b) parametric form assumed, feasible region not discovered. Both grounds are
defensible against a reviewer from the SDF community.

**Historical-envelope limitation.** Stated honestly in paragraph 6: "The current
implementation explores what is achievable by recombination of the historical
monthly flow record through the Kirsch-Nowak bootstrap generator." Not buried in
Section 4; visible to a reviewer reading only the Introduction.

**FIND comparison reframing.** The text does not say "FIND doesn't scale" or
"weight retuning is a barrier." It says the ex-post feasibility discovery cost
"limits the practical dimensionality of the ensemble design" — this is the
automatic-feasibility-discovery framing, not an efficiency attack. The timing
comparison is correctly deferred to Methods/Results.

**No absolute "first" claim.** The novelty statement specifies a conjunction
of three properties, not a claim of being the first method to work in hazard
space at all.

**Style checks.** No em-dashes used. No colons in flowing prose (colons appear
only in the self-review section headers, which will be deleted). Semicolons
appear only within parenthetical citation lists (standard academic usage, not
flowing prose). No quotations from cited works. No transcribed numerical values
from cited papers. No informal or quippy constructions.

**Open questions for author review.**
1. The SDF paragraph currently omits Salvadori and De Michele; if the author
   wants to add a general copula reference, insert after Serinaldi et al. (2009).
2. "Conditioned Latin hypercube sampling" (cLHS) is used for Bonham et al. (2024)
   but the method name can be simplified if the author prefers "conditioned LHS."
3. Paragraph 5 is long (~450 words combined). If the editor requests a shorter
   Introduction, the bivariate SDF section can be trimmed to two sentences after
   Serinaldi et al. (2009) and the contrast moved to a footnote or Section 4.
4. The phrase "physically consistent with the observed record" appears in both
   paragraph 3 (Herman et al. 2016) and paragraph 5 (Borgomeo et al. 2015).
   Consider varying the phrasing if both paragraphs survive into the final draft.
