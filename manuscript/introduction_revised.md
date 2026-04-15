# MOEA-FIND Introduction — Revised Draft (Round 3, post-critique-2)

*Branch: claude/zen-swanson. Two critique/revise cycles complete (autonomous session 2026-04-14).*
*Round 1 backup: introduction_revised_round1.md.*
*Round 2 backup: introduction_revised_round2.md.*
*Critiques: intro_critique_round1.md, intro_critique_round2.md.*
*Key open items for Trev listed in self-review notes below.*

---

## 1. Introduction

Multi-year meteorological drought is the binding supply-side risk for reservoir-based
water supply systems in water-stressed regions, and robust evaluation of operating
policies and drought management rules under deep uncertainty requires synthetic
streamflow ensembles that span a meaningful range of physically plausible drought
conditions. The many-objective robust decision making framework of Kasprzyk et al.
(2013), building on the robust decision making tradition of Lempert et al. (2003)
and the many-objective vulnerability analysis of Herman et al. (2015), evaluates
candidate operating policies by simulating each across an ensemble of plausible
future states of the world and aggregating system performance into robustness
measures that guide policy selection. The conclusions a planner can draw about
policy vulnerability depend directly on how thoroughly that ensemble spans the
drought hazard space defined by the severity, duration, seasonal timing, and
system-specific threshold characteristics that govern the operating rules and
robustness metrics for the basin under study. The exploratory modelling, synthetic
hydrology, and scenario discovery communities have developed a portfolio of
complementary approaches to synthetic scenario design, each concentrating the
structured design effort at a different stage of the analysis pipeline, as shown
in Figure 1. Bonham et al. (2025) demonstrate that ensembles designed for coverage
in the uncertain input space systematically undersample challenging conditions in the
downstream performance space, confirming that input-space design criteria alone do
not guarantee adequate coverage of the drought hazard conditions that matter most
for policy evaluation. This paper introduces a fourth approach, MOEA-FIND, that
targets the drought hazard outcome space directly through multi-objective
evolutionary search rather than as an emergent projection of a design applied
upstream of that space.

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
specifying them in advance.

A third approach, library-and-subsample, generates a large pre-computed library of
synthetic traces from an unsteered stochastic generator and applies a structured
subsampling criterion to select a tractable ensemble for downstream simulation
(Figure 1, panel C). This approach shares the forward-propagation logic of panel A
but makes a qualitatively distinct methodological choice at the second stage.
Rather than embedding the coverage criterion implicitly in the input sampling
distribution before any synthetic traces are generated, library-and-subsample makes
the coverage criterion an explicit, inspectable decision applied within the
coordinates of the pre-generated library, so the relationship between the design
choice and the resulting ensemble is visible and adjustable after generation. Bonham
et al. (2024) developed this approach for the Colorado River Basin by generating a
multi-scenario library of streamflow, demand, and initial storage combinations and
applying conditioned Latin hypercube sampling to select a subsample that covers the
generator input coordinates with near-uniform spacing. Their evaluation demonstrates
that input-space space-filling metrics are strong predictors of robustness ranking
preservation as subsample size decreases, so that a carefully selected small
ensemble recovers the ranking structure of the full library when the subsampling
criterion is well matched to the input space. This approach is computationally
attractive because the generator is run once to produce the library and the
subsampling step is inexpensive relative to regenerating traces or running the
downstream simulation model. The limitation that applies to library-and-subsample
methods is that the subsampling criterion is defined in the input coordinates of the
generator rather than in the drought hazard characteristic space the vulnerability
analysis evaluates. The generator map from input coordinates to drought
characteristics is not in general monotone or invertible, so a space-filling design
in the input space is expected to produce a non-uniform projection into the hazard
space, oversampling the drought conditions the generator produces most readily from
the historical record and undersampling the tails and combinations that are
physically attainable but statistically rare. A quantitative characterization of
this non-uniformity for the Kirsch-Nowak generator applied to the Cannonsville basin
is provided in Section 3. Bonham et al. (2025) confirm the projection distortion
as a systematic bias in input-space ensembles across multiple case studies. The
analyst who wants structured coverage of the drought hazard space must either accept
the emergent projection of the input-space design or apply a secondary subsampling
step in the hazard coordinates, which requires a library large enough that the
hazard space is already well populated before subsampling.

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
methods for joint extreme event frequency analysis in hydrology that has been widely
applied to the joint distribution of drought event attributes in subsequent work.
Copula-based severity-duration-frequency methods are valid and well-validated for
probabilistic drought hazard characterization, and they provide a principled way to
generate scenarios with specified joint return periods of drought severity and
duration. The limitation relevant to the present paper applies to the parametric-
form assumption that is shared by all such methods, whether operating as standalone
event-attribute samplers or as components within a conditional streamflow
disaggregation model. In all copula-based severity-duration-frequency formulations,
a parametric joint distribution governs which drought signatures are generated, and
that distribution is fitted to the historical event record rather than discovered
from the generator by directed search. MOEA-FIND does not assume a parametric form
for the joint distribution of drought characteristics. It discovers the achievable
region from the generator empirically, so the shape and extent of the feasible
drought hazard region is a byproduct of the optimization rather than a modelling
assumption imposed before the search begins.

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
with secondary terms penalizing non-drought-period flow statistics including mean
annual flow and flow percentiles. MOEA-FIND externalizes equivalent non-drought
plausibility controls into a separate constraint layer rather than objective terms,
which preserves the orthogonality of the drought characteristic objectives and
allows the Manhattan-distance construction to function. The design choice is made
explicit in Section 2. Wheeler et al. (in press) extended the single-objective
target-matching formulation to multi-site basins, expanding the target vector with
cross-site correlation matrices and producing one joint multi-site record per
high-dimensional target vector. The structural limitation common to all
target-matching methods is that each optimization run produces one synthetic trace
that matches one pre-specified target vector, so building an ensemble that spans a
range of drought characteristic combinations requires running the optimizer once per
target specification. The subset of the target combination space that is physically
attainable under the generator is not known in advance, and an analyst enumerating
a target grid must verify attainability for each cell independently. At the three-
characteristic scale demonstrated by Zaniolo et al. (2024) this step is tractable,
but extension to higher-dimensional characteristic spaces would require verifying a
grid that grows exponentially in the number of cells as characteristics are added.
The gap that motivates the present paper is the absence of a method that both
discovers the feasible drought hazard region automatically from an arbitrary
stochastic generator and delivers structured coverage of that region in a single
optimization run, without requiring the analyst to enumerate target combinations in
advance and without assuming a parametric form for the joint distribution of
achievable drought characteristics.

This paper introduces MOEA-FIND, a fourth approach that places the structured
design step directly in the drought hazard characteristic space (Figure 1, panel D).
The method explores which drought signatures are achievable by recombination of
observed monthly values through the Kirsch-Nowak bootstrap generator (Kirsch et al.,
2013; Nowak et al., 2010) and does not extrapolate beyond the historical flow
envelope. Parametric extension for climate-change stress testing, through
substitution of a kappa-distribution generator or vine-copula model at the same
coupling point, is architecturally compatible and is deferred to follow-on work.
Within that scope, MOEA-FIND couples the Borg multi-objective evolutionary algorithm
(Hadka and Reed, 2013) to the Kirsch-Nowak generator through a K+1-objective
formulation in which a single auxiliary Manhattan-distance objective causes every
feasible synthetic trace to be Pareto-optimal, not merely those on the
non-dominated frontier of the drought characteristic objectives alone. This
formulation converts coverage of the feasible drought hazard region from a post-hoc
subsampling step into a property that follows from the structure of the epsilon-
dominance archive, given sufficient optimizer convergence. The archive tiles the
feasible region at a resolution determined by the user-supplied epsilon vector, and
convergence of the archive constitutes the structured ensemble that is the
deliverable of the method. Two design properties require explicit statement. First,
MOEA-FIND does not prescribe what drought characteristics a basin can produce and
does not redefine what a drought is. The Kirsch-Nowak generator and its historical
forcing remain the physical and statistical model of drought emergence. The optimizer
determines which combinations of drought characteristics are achievable under that
model without imposing a parametric form on their joint distribution. Second,
physical correlations among drought characteristics, such as the tendency for longer
events to accumulate larger cumulative deficits and therefore to be more severe on
average, are a feature of the feasible drought hazard region rather than a
limitation of the method. The method is designed so that the Pareto-optimal archive
adapts to the correlated structure of the feasible region and covers its interior,
including portions where physically coupled characteristics co-occur. Section 3
evaluates empirically whether this coverage is achieved on the Kirsch-Nowak
generator applied to the Cannonsville basin. Marginal and pairwise diagnostics for
the archive members are reported in the coverage analysis as a check on this
property. The novelty of the present contribution is the K+1-objective formulation
coupled with an epsilon-dominance archive that tiles the resulting hyperplane at a
user-specified resolution. No existing method uses an auxiliary Manhattan-distance
objective that causes all feasible synthetic traces to be Pareto-optimal by
construction combined with an epsilon-dominance archive that converts this property
into structured coverage of the feasible hazard region in a single optimization run.
Prior approaches to hazard-space targeting either fit parametric distributions to the
event record and sample from the fitted model, as in the copula-based severity-
duration-frequency tradition, or optimize one synthetic trace per user-specified
target and discover the feasible region ex post, as in the FIND tradition. MOEA-FIND
combines the discovery and coverage steps in a single optimization run that requires
no prior specification of attainable target combinations and no parametric assumption
about the joint distribution of the drought characteristics.

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
generator map is known, with the analog verification on the Kirsch-Nowak application
reported alongside multi-seed convergence diagnostics and library-subsample coverage
comparisons in Section 3. The fourth is a demonstration of the method on the
Cannonsville inflow with three custom drought characteristics, including a cyclic
peak timing feature, alongside a coverage comparison against equal-size library-
subsample baselines and a scenario discovery exercise conducted in drought hazard
space using operational failure labels derived from Pywr-DRB reservoir simulation
as the classification target. Section 2 defines the drought characterization
pipeline, the structured ensemble generation method, the analytic benchmark problem,
the Cannonsville case study, and the coverage and scenario discovery evaluation
protocols. Section 3 presents the three sets of empirical results. Section 4
discusses limitations, scope, and extensions. Section 5 concludes.

---

*End of introduction body. Self-review notes follow — delete before submission.*

---

## Self-review notes — round 2 critique responses

**Issue 1 (epsilon-NSGA-II loophole — carry-over):** Addressed. Novelty sentence
now reads "No existing method uses an auxiliary Manhattan-distance objective that
causes all feasible synthetic traces to be Pareto-optimal by construction combined
with an epsilon-dominance archive that converts this property into structured
coverage of the feasible hazard region in a single optimization run." The
formulation-plus-archive conjunction is now named explicitly. The claim is that no
prior PUBLICATION has used this combination, not that no other algorithm could be
configured to do so.

**Issue 2 (hybrid copula regression):** Addressed without conceding gap scope.
Rather than saying "hybrids do exist and do not share this limitation" (which
narrowed the gap), the revised text says the parametric-form assumption is shared
"by all such methods, whether operating as standalone event-attribute samplers or
as components within a conditional streamflow disaggregation model." This is a
defensible universal claim: all copula-based formulations impose a parametric joint
distribution, regardless of coupling. The "hybrid" literature concession is removed.
If a reviewer from the hybrid SDF community objects, the response is: even in the
coupled case, the parametric copula governs which drought signatures are generated,
and this assumption is what MOEA-FIND eliminates. A citation for the hybrid class
was not added because no specific verifiable paper could be confirmed; if Trev wants
to cite one, it should go here.

**Issue 3 (library-and-subsample separation — carry-over):** Addressed with a
stronger epistemic argument. The new sentence reads: "Rather than embedding the
coverage criterion implicitly in the input sampling distribution before any synthetic
traces are generated, library-and-subsample makes the coverage criterion an explicit,
inspectable decision applied within the coordinates of the pre-generated library, so
the relationship between the design choice and the resulting ensemble is visible and
adjustable after generation." This argues an epistemic difference (implicit pre-
generation vs. explicit post-generation visibility) rather than only a procedural
one (before vs. after). Whether this convinces a hostile reviewer is debatable, but
the argument is now stated rather than assumed.

**Issue 4 (Bonham placement regression):** Addressed. Bonham et al. (2025) is now
at sentence 5 of Para 1, positioned as the bridge between the four-approaches
taxonomy (sentence 4) and the MOEA-FIND introduction (sentence 6). The logical
flow is: motivation (sentences 1-3) → existing portfolio of approaches (sentence 4)
→ empirical evidence that the portfolio leaves a gap (sentence 5) → this paper fills
it (sentence 6). The flow is intact.

**Issue 5 (guaranteed — carry-over):** Addressed. "Guaranteed by the structure of
the epsilon-dominance archive" replaced with "follows from the structure of the
epsilon-dominance archive, given sufficient optimizer convergence." This preserves
the architectural claim while acknowledging that convergence to the interior is what
Section 3 evaluates empirically.

**Issue 6 (FIND convergence-failure claim — new):** Addressed. The claim about
convergence failure as an infeasibility-detection mechanism is removed. The revised
text says "the subset of the target combination space that is physically attainable
under the generator is not known in advance, and an analyst enumerating a target
grid must verify attainability for each cell independently. At the three-
characteristic scale demonstrated by Zaniolo et al. (2024) this step is tractable,
but extension to higher-dimensional characteristic spaces would require verifying a
grid that grows exponentially in the number of cells as characteristics are added."
This characterizes the limitation as prospective scaling cost rather than as a
convergence failure that FIND already exhibits.

**Issue 7 (exponential grid as prospective — new):** Addressed jointly with Issue 6.
The exponential-grid argument is now framed as a prospective concern, not a
demonstrated limitation.

**Issue 8 (style violations — new):** All six identified style violations have been
corrected. The text has been scanned for colons in flowing prose (none remaining)
and semicolons in flowing prose (none remaining other than within parenthetical
citation lists, which are permitted). Complete check:
- Colon at "standalone severity-duration-frequency models: these methods..." → Period
  replaced; text restructured as two sentences.
- Semicolon at "do not share this limitation; the distinction..." → Removed by
  restructuring the hybrid-copula handling (Issue 2).
- Semicolon at "construction to function; the design choice..." → Split into two
  sentences.
- Semicolon at "model of drought emergence; the optimizer determines..." → Split
  into two sentences.
- Semicolon at "characteristics co-occur; Section 3 evaluates..." → Split into two
  sentences.
- Semicolon at "generator map is known; the analog verification..." → Rephrased as
  "with the analog verification... reported alongside..." in a relative clause.

**Issue 9 (non-monotone assertion — residual):** Addressed. Added sentence: "A
quantitative characterization of this non-uniformity for the Kirsch-Nowak generator
applied to the Cannonsville basin is provided in Section 3." This anchors the
logical claim empirically without asserting numbers that belong in Section 3.

**Remaining open items for Trev (in priority order):**

1. **Wheeler et al. (in press):** Verify complete bibliographic entry, venue, and
   DOI or preprint URL. If not yet in review, either remove the citation or replace
   with a factual description of the multi-site extension without attribution.

2. **Bonham et al. (2025):** Verify full title and DOI (distinct from 2024 paper).
   Both PDFs are in the literature/ directory, but the 2025 paper's exact title and
   citation should be confirmed from the source.

3. **Salvadori and De Michele (2004):** Confirm this is the correct year and venue
   (expected: Water Resources Research, 2004, vol. 40, W12511 or similar). The
   citation is used in the SDF paragraph without a page or doi check.

4. **Hybrid SDF-generator citation:** If Trev wants to acknowledge the hybrid
   literature with a specific citation rather than the general claim, add it after
   Salvadori and De Michele (2004). The current text makes a defensible universal
   claim without the citation; adding one would strengthen it.

5. **Fourth contribution (Pywr-DRB labels):** The Section 3.3 redesign (operational
   failure labels from Pywr-DRB simulation) requires HPC Phase gamma runs to be
   complete before submission. If those runs are not complete, revise the fourth
   contribution to frame it as a workflow demonstration with proxy failure labels,
   with operational validation deferred to the companion paper.

6. **"Fourth approach" taxonomy:** The round-2 critique (Issue 3) notes this
   remains "adequately defended but still vulnerable." A one-paragraph review
   response is implicit in the Para 4 text, but Trev should be prepared to defend
   it explicitly if a reviewer requests collapsing panels A and C.
