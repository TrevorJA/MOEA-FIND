# MOEA-FIND: A Multi-Objective Optimization Framework for Structured Synthetic Drought Ensemble Generation

*SCAFFOLD DRAFT. Section headers and a small amount of literature-grounded
motivation only. Every other section is intentionally a placeholder until the
HPC empirical work has been completed and the corresponding decision in
`governance/design_decisions.md` is settled. Do not draft prose into a section
whose underlying result is not yet in hand.*

*Working draft. Target venue: Water Resources Research. Target length: approximately 10,500 words main text, seven figures and one table, separate Supporting Information.*

Authors (provisional): Trevor Amestoy and co-authors to be determined.

Last updated: 2026-04-27 (scaffold reset).

---

## Key Points

> *Placeholder. Three Key Points to be drafted after Sections 3.1, 3.2, and
> 3.3 contain final empirical results and after the contribution framing has
> been validated against the production Cannonsville archive.*

## Plain Language Summary

> *Placeholder. To be drafted last, after the Abstract is finalised.*

## Abstract

> *Placeholder. The Abstract describes results that have not yet been
> produced and must wait until the empirical sections are populated.
> Premature abstract prose is the most likely source of overclaim and is
> avoided here by design.*

---

## 1. Introduction

Water supply planning unfolds across decades during which climate, demand,
and institutional conditions can change in ways that cannot be represented by
any single probabilistic forecast. Methods for evaluating reservoir operating
policies and drought management rules under this kind of deep uncertainty
(Walker et al., 2013) have converged on robustness-based frameworks in which
a candidate policy is judged not by its expected performance under a
best-estimate projection but by how its performance varies across an ensemble
of plausible future conditions (Lempert et al., 2003; Kasprzyk et al., 2013;
Herman et al., 2015). Moallemi et al. (2020b) characterise the computational
workflow underlying these frameworks as a mapping from a decision space of
candidate alternatives and an uncertainty space of plausible states of the
world into an outcome space of system performance measures, where robustness
is inferred from the distribution of outcomes across the sampled uncertainty
space (Moallemi et al., 2020a). The quality of any robustness inference drawn
from this workflow is bounded by how well the uncertainty sample represents
the conditions that drive policy performance for the system under study. For
water supply systems in which multi-year meteorological drought is a
principal source of operational stress, the uncertainty sample takes the form
of an ensemble of streamflow records, and the conditions that drive policy
performance are most usefully described in terms of the drought
characteristics the ensemble contains.

Ensembles of plausible streamflow records for robustness analyses are
constructed using synthetic streamflow generators, which draw stochastic
realisations of monthly or daily flow from a model calibrated to the
historical record. The design target shared across generator families is to
produce multiple traces whose statistical properties are indistinguishable
from those of the observed record on the features the downstream analysis
relies upon. Parametric autoregressive models of the Thomas-Fiering type
(Thomas and Fiering, 1962) and their seasonal extensions impose these
properties through an explicit distributional assumption, while nonparametric
approaches including moving-block bootstrap and nearest-neighbour resampling
instead preserve them through repeated sampling of the observed record.
Kirsch et al. (2013) introduced a modified fractional Gaussian noise procedure
that performs an uncorrelated bootstrap resampling of historical monthly flows
followed by Cholesky correlation restoration to impose the historical monthly
correlation structure in Gaussian space, an approach that has been adopted
widely in the water supply vulnerability analysis literature (Hadjimichael et
al., 2020a, 2020b; Quinn et al., 2018). Herman et al. (2016) developed a
weighted bootstrap extension that emphasises the low-flow quantiles of the
historical record and yields ensembles with analyst-controlled increases in
drought frequency. Parametric alternatives based on kappa-distribution
marginals and copula dependence structures can extend the synthesised drought
envelope beyond what historical recombination supports at the cost of
committing to an explicit distributional form (Svensson et al., 2017;
Brechmann et al., 2017). A related tradition in drought frequency analysis
characterises the joint probability distribution of drought severity and
duration through bivariate or multivariate copula models fitted to historical
or regional event records (Shiau, 2006; Serinaldi et al., 2009); these
severity-duration-frequency methods are designed for return-period estimation
and risk quantification rather than for the generation of synthetic streamflow
traces that can be propagated through a water system simulation model.

The dominant practice for constructing a robustness ensemble samples deeply
uncertain factors from a space of plausible states of the world using a
space-filling experimental design such as Latin hypercube sampling (McPhail
et al., 2018; Hadjimichael et al., 2020a; Quinn et al., 2020). The sampled
factors typically include hydroclimatic perturbations applied to the
historical record, hyperparameters of the streamflow generator itself, and
human-system drivers such as demand growth, infrastructure commitments, and
institutional constraints. Following the nested-loop structure formalised by
Hadjimichael et al. (2020a), each sampled state of the world indexes an
outer condition under which an ensemble of inner stochastic realisations is
drawn from the generator, and the union of outcomes across all realisations
is interrogated through scenario discovery (Bryant and Lempert, 2010;
Hadjimichael et al., 2020a), sensitivity analysis, or direct robustness
scoring. The space-filling design ensures approximately uniform coverage of
the upstream uncertainty factors, and forward propagation through the
generator and simulation model produces the outcome ensemble from which
robustness is inferred.

This workflow has been productive for identifying which uncertain factors
condition system vulnerability, but it has three connected limitations that
bear directly on drought-driven vulnerability analysis. The first is that
the mapping from the upstream uncertainty space through the generator and
the simulation model into the drought characteristic space that drives
operational outcomes is nonlinear and not one-to-one, so that a sample that
is approximately uniform in the upstream space need not produce approximately
uniform coverage of the drought characteristic space that governs policy
performance. Quinn et al. (2020) demonstrate that the vulnerability
conclusions drawn from exploratory modelling depend on the analyst's choice
of sampling design and assumed probability distributions, and argue that
this sensitivity should itself be assessed. The second limitation is that
drought combinations causing operational failure typically occupy a small
fraction of the drought characteristic space relative to benign and
historically common conditions, so that under realistic computational budgets
the rare but operationally consequential combinations are systematically
under-represented while common and non-stressful conditions are
over-represented. Bonham et al. (2025) characterise this as a general
property of ensembles derived from climate projections and historical
resampling, and note that it produces complex and difficult-to-anticipate
impacts on the factor-mapping steps that follow. The third limitation is
that much of the simulation effort in an upstream-sampled ensemble is spent
on traces during which no operationally relevant drought occurs, yielding
simulation records that are uninformative for the vulnerability assessment
the analysis is designed to produce. Bonham et al. (2024) ground the severity
of this computational constraint in the context of a Colorado River Basin
robustness analysis and propose a post-hoc statistical subsampling of the
precomputed ensemble that preserves policy rankings at reduced cost, but
their subsampling operates on the same uncertainty factor axes that
originally defined the ensemble and does not alter the relationship between
the upstream sampling design and the drought characteristic space.

A separate line of work addresses the parameter-to-outcome mapping problem
from the opposite direction, coupling a streamflow generator directly to an
optimisation procedure that searches for synthetic traces matching
analyst-specified drought characteristics. Borgomeo et al. (2015) introduced
a simulated annealing procedure in which the decision vector is a permutation
of bootstrapped monthly historical flows and the objective function is a
weighted sum of deviations from a user-specified vector of hydrological
targets, producing one synthetic trace per optimisation run that matches the
supplied target while remaining consistent with the observed record.
Borgomeo et al. identified the extension of the formulation to multi-objective
search as an open direction. Zaniolo et al. (2024) adapted the Borgomeo
procedure specifically to drought targeting through the FIND method, which
defines three drought characteristics on the Standardized Streamflow Index
(frequency, intensity, and duration) and minimises a weighted-sum deviation
from user-specified values of each characteristic; the FIND method requires
per-regime weight retuning when targets move away from historical conditions,
and Zaniolo et al. demonstrate the approach on a five-by-five grid of
intensity and duration increments, each requiring an independent optimisation
run. Wheeler et al. (2025) extended the single-objective target-matching
formulation to multi-site basins by expanding the target vector with cross-
site correlation matrices and a Hurst coefficient for long-range dependence.
These contributions share a common and productive shift in framing relative
to the upstream-sampling tradition, moving the analyst's specification of
drought conditions from the uncertainty factor space of the generator to the
drought characteristic space that the vulnerability analysis evaluates. They
also share a common structural constraint: each optimisation run produces
one synthetic trace corresponding to one pre-specified target vector, and an
ensemble spanning a range of drought conditions is assembled by running the
optimiser repeatedly against an analyst-enumerated grid of target vectors.
The physically attainable portion of that grid is not mapped before search,
and the enumeration cost grows rapidly with the number of drought
characteristics the analysis requires.

The upstream-sampling tradition and the search-based tradition therefore
face complementary limitations for drought ensemble construction. The former
produces structured coverage of the uncertainty factor space but leaves
coverage of the drought characteristic space as an emergent and uncontrolled
property of the forward mapping. The latter operates directly in the drought
characteristic space but requires the analyst to enumerate the targets the
search will match, without advance knowledge of which targets are attainable.
The design problem this paper addresses is how to produce, under a realistic
computational budget, an ensemble of synthetic streamflow traces whose
coverage is structured directly in the drought characteristic space that the
vulnerability analysis evaluates and whose member traces recover the feasible
drought characteristic region as a byproduct of the generation step.

> *Remainder of §1 (statement of contribution, method overview, paper
> roadmap) deferred until the contribution can be supported by Sections 3.1,
> 3.2, and 3.3.*

---

## 2. Methods

> *Section 2 is a structural placeholder. The Manhattan-distance auxiliary
> objective construction in §2.2.3 is the only methodological commitment at
> this stage; every other subsection (drought characterisation, generator
> details, optimiser settings, benchmark problem, case-study setup, coverage
> and evaluation protocols) is under active investigation and will be
> drafted only when the supporting decisions in
> `governance/design_decisions.md` and the corresponding HPC results are
> settled.*

### 2.1 Quantification of drought characteristics

> *Pending prose. The production metric set is now settled
> (DD-04, 2026-04-27): mean event severity, mean event cumulative
> deficit, and time-in-drought fraction — three continuous metrics
> spanning depth, volume, and persistence axes. Drought duration and
> peak-severity month are excluded for clustering at discrete monthly
> values. Reference: `governance/design_decisions.md` §DD-04;
> implementation in `src/drought_metrics.py` (`primary` preset) and
> `src/objectives.py`. Detailed prose deferred until §3 results
> motivate the framing.*

### 2.2 Generation of structured drought hazard ensembles

#### 2.2.1 Kirsch-Nowak synthetic streamflow generator

> *Pending. References: Kirsch et al. (2013); Nowak et al. (2010); Herman et
> al. (2016); `src/kirsch_wrapper.py`.*

#### 2.2.2 Borg MOEA and epsilon-dominance archiving

> *Pending. References: Hadka and Reed (2013); Reed et al. (2013); Laumanns
> et al. (2002); `governance/design_decisions.md` §DD-07.*

#### 2.2.3 Manhattan-distance auxiliary objective

> *Pending. The construction is the load-bearing methodological contribution
> of the paper. The L1 form ($f_j = D_j$, $f_{K+1} = \lVert D - D^* \rVert_1$)
> is fixed and implemented in `src/objectives.py`; see
> `governance/design_decisions.md` §DD-11. The full derivation will be
> drafted into Supporting Information once the empirical sections that
> motivate the framing have been written.*

#### 2.2.4 Algorithm and implementation

> *Pending. Reference: `planning/code_state.md` for current pipeline.
> Specific NFE budget, epsilon vector, and plausibility constraint regime
> are pending HPC results and the resolution of `code_alignment_backlog.md`
> Items 3, 8.*

### 2.3 Analytic benchmark problem

> *Pending. The choice of analytic benchmark, dimensionalities tested, NFE
> budget, and seed protocol are not yet final. Preliminary dimension-sweep
> evidence on a constrained K-dimensional benchmark is in
> `evidence/shell_vs_interior_diagnostic.md` and DD-11.*

### 2.4 Delaware River Basin case study

> *Pending. The Cannonsville single-site setup is the case study for this
> paper. Multi-site DRB extension is deferred to a follow-up paper per
> `planning/publication_plan.md` §4. Anti-ideal placement protocol pending
> the verification described in DD-11 anti-ideal assumption.*

### 2.5 Coverage metrics and evaluation protocols

> *Pending. The Section 3.3 failure label has been redesigned to use a
> Pywr-DRB-derived operational outcome rather than a threshold on the same
> drought characteristics that MOEA-FIND optimises over (HC-2 resolution).
> The full protocol is in `planning/experiment_plan.md` Part A.*

---

## 3. Results

> *Section 3 is a placeholder. No prose belongs here until the corresponding
> HPC results are in hand.*

### 3.1 Interior coverage on the analytic benchmark problem

> *Pending. Numerical results in
> `evidence/shell_vs_interior_diagnostic.md` are produced by MM Borg MOEA
> on HPC (per DD-07, updated 2026-04-28). Drafting blocked only on the
> final figure pass.*

### 3.2 Coverage of the feasible drought hazard region at Cannonsville

> *Pending HPC Phase β (single-site Cannonsville) and Phase γ (10 k-trace
> library baseline). Drafting is gated on DD-12.*

### 3.3 Scenario discovery in the drought hazard space

> *Pending HPC Phase β and γ. The failure label will be a Pywr-DRB-derived
> operational outcome (HC-2 resolution); see `planning/experiment_plan.md`
> Part A. Specific drought-characteristic features used as classifier inputs
> are not yet final.*

---

## 4. Discussion

> *Pending. Discussion prose depends on the empirical results in §3 and on
> the production Cannonsville archive. Drafting any of the limitations,
> scope, or extension paragraphs before the results land risks anchoring the
> discussion to outcomes that the empirical work has not yet supported.*

---

## 5. Conclusions

> *Pending. To be drafted last, in approximately 400 words, after §§3 and 4
> have been written.*

---

## Open Research

> *Placeholder. Code and data availability statement to be finalised at
> submission. Current code state is summarised in `planning/code_state.md`.*

## Acknowledgments

> *Placeholder. HPC resources, co-authors, and funding sources to be listed
> at submission.*

## References

> *Placeholder. Bibliography to be assembled at submission. The
> literature-anchored Introduction above cites Walker et al. (2013); Lempert
> et al. (2003); Kasprzyk et al. (2013); Herman et al. (2015); Moallemi et
> al. (2020a, 2020b); Thomas and Fiering (1962); Kirsch et al. (2013); Nowak
> et al. (2010); Herman et al. (2016); Hadjimichael et al. (2020a, 2020b);
> Quinn et al. (2018, 2020); McPhail et al. (2018); Bryant and Lempert
> (2010); Bonham et al. (2024, 2025); Borgomeo et al. (2015); Zaniolo et al.
> (2024); Wheeler et al. (2025); Svensson et al. (2017); Brechmann et al.
> (2017); Shiau (2006); Serinaldi et al. (2009).*
