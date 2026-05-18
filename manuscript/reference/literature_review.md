# Comprehensive Literature Review: MOEA-FIND

*Date: 2026-04-12. Updated 2026-04-14 after the opus predecessor-PDF
critique corrected stale Zaniolo and Wheeler citations, confirmed the
target-matching characterisation of all three predecessors, and
corrected the Wheeler journal reference from Journal of Water
Resources Planning and Management to Journal of Hydrologic Engineering.
The K = 2 through K = 6 shell-versus-interior dimension sweep
documented in `evidence/shell_vs_interior_diagnostic.md` and in DD-11 is
the empirical anchor for the corrected method exposition and is not
reflected in the pre-2026-04-14 entries below.*
*Purpose: Ground all experimental design and method choices in published literature*

---

## 1. Search-Based Synthetic Streamflow Generation

### Borgomeo et al. (2015) — "Numerical Rivers"
**Citation:** Borgomeo, E., et al. (2015). Numerical rivers: A synthetic streamflow generator for water resources vulnerability assessments. Water Resources Research, 51(7).

Pioneered treating streamflow generation as a combinatorial optimization problem. Uses simulated annealing to shuffle values in initialized synthetic time series until sequences match user-defined objectives (target statistical properties). Enables representation of climate-induced changes in persistence, variability, and extremes. Single-objective formulation: one target, one trace per optimization run.

**Relevance:** Direct intellectual predecessor. MOEA-FIND generalizes to multi-objective, where the Pareto front is the ensemble.

### Zaniolo et al. (2024) — FIND
**Citation:** Zaniolo, M., Fletcher, S., and Mauter, M. (2024). FIND: A synthetic weather generator to control drought Frequency, INtensity, and Duration. Environmental Modelling and Software, 172, 105927.

Extends Borgomeo's approach with direct control of three drought properties (frequency, intensity, duration) via the Standardized Streamflow Index, along with autocorrelation and non-drought percentile penalty terms, all aggregated into a single weighted-sum simulated annealing objective. The user supplies a single `(F, I, D)` target per run, and the block-replace-from-distribution operator draws replacement values from historical `n_months`-cumulative CDFs and disaggregates via k-NN. The Experiment 3 demonstration enumerates a `5 x 5` grid of intensity and duration increments and runs the algorithm 25 times, three replicates per cell. Per-regime weight retuning is explicitly required when targets move away from historical conditions.

**Relevance:** Closest precedent. MOEA-FIND replaces FIND's single-objective target matching with a multi-objective feasibility discovery formulation that returns the grid of attainable targets as the deliverable of one run, so no per-target weight retuning is necessary.

### Wheeler et al. (2025) — Multi-Site Extension
**Citation:** Wheeler, K.G., Simpson, M., Borgomeo, E., and Hall, J.W. (2025). Multisite Nonparametric Stochastic Streamflow Generation for the Eastern Nile Basin. Journal of Hydrologic Engineering, 30(1), 04024056.

Extends the Borgomeo (2015) simulated annealing algorithm to multiple sites by expanding the target vector with cross-site correlation matrices and swapping values at a single randomly selected site per iteration. Single-objective weighted-sum formulation is preserved. The target vector includes mean monthly flow, monthly standard deviation, monthly autocorrelation to lag 11, interannual standard deviation, interannual lag-1 autocorrelation, the Hurst coefficient, and two cross-correlation matrices at each of 18 gauges. The 100-trace ensemble is 100 replicate runs against one historical target vector. The climate-change demonstration uses a single Hurst perturbation factor `p = 1.2` and leaves exploration of broader perturbations to future work.

**Relevance:** Demonstrates that search-based target matching scales to multi-site basins and demonstrates Hurst coefficient inclusion as a long-range dependence constraint. MOEA-FIND uses shared bootstrap indices across sites as a simpler cross-site correlation mechanism and does not place cross-site structure in the optimization objective.

### Gozini et al. (2026) — Parametric Inverse Approach for Multi-Site Stress Testing
**Citation:** Gozini, H., Asadzadeh, M., Stadnyk, T., and Slota, P. (2026). Inflow generation for water system planning in multisite studies using target synthetic streamflow generation. Water Resources Research, 62, e2025WR042355.

Closes the gap left by Borgomeo (2015) and Wheeler (2025) by applying the inverse approach to a *parametric* generator (Kirsch et al. 2013) instead of a non-parametric one. The decision variables are forcing percent changes `(M_F, SD_F)` per (month, location); the optimization is decoupled into one 2-D scalar L1 minimization per (month, location) cell using SciPy `differential_evolution`. The target tensor is constructed multiplicatively from a per-(month, location) seasonality layer and a per-scenario uniform deviation layer (Eq 1). The 14-location Winnipeg River Basin demonstration sweeps a 17 × 17 grid of deviations from `-80%` to `+80%` in `10%` increments (289 scenarios). A closed-form lognormal feasibility polygon, computed per (month, location), is enforced upstream of the optimization by snapping infeasible target cells onto the polygon boundary. Cross-site correlation is preserved by sharing the matrix year across all locations within a scenario, identical to the Quinn et al. (2018) Kirsch-Nowak convention. Daily disaggregation uses Nowak et al. (2010) k-nearest neighbors. Reported runtime is 3.22 hr / 100 scenarios on a 12-core laptop, compared to 169 hr for Wheeler (2025) and 4 hr for the single-site Borgomeo (2015).

The paper provides the most thorough statistical-credibility audit in this lineage: Wilcoxon rank-sum and Levene's tests on monthly distributions, autocorrelation function to 50 months lag with RMSE against historical, pairwise cross-correlation RMSE, and the same battery applied separately to a stationary case (RH-O), an expansion-over-historical case (EH-O), and a GCM-anchored seasonality-shift case (EG-O). The stress-testing demonstration runs three performance metrics (annual + seasonal hydropower reliability, flood frequency, water-demand shortage) on the MB_HydroSim Manitoba Hydro operations model and produces success/failure response surfaces in 2-D (mean × SD) deviation space.

**Relevance — common ground:** Same Kirsch (2013) generator core. Same Nowak (2010) k-NN disaggregator. Same shared-matrix-year cross-site correlation mechanism. Same paradigm-(iv) directed-search lineage as Borgomeo, Zaniolo, Wheeler. Same six-requirement framing for plausible scenarios (statistical credibility, target-property control, deep-uncertainty range, uniform exposure-space coverage, physical and mathematical feasibility, computational efficiency) that MOEA-FIND inherits from this lineage. Same explicit feasibility framing (their Figure 3 versus MOEA-FIND DD-10).

**Relevance — distinction:** Gozini's exposure space is the 2-D plane of monthly mean and SD deviations per (month, location) — first/second moments of the marginal distribution at each cell. Their paper explicitly grounds this in Quinn et al. (2018) Figure 9, which they cite as the canonical hydrologic-stressor exposure space (their §1, page 2). MOEA-FIND's exposure space is the K-dimensional emergent drought-hazard characteristic space (SSI event severity, cumulative deficit, time-in-drought fraction; DD-04). The distinction is not cosmetic: their per-cell decoupled optimization works *only because the targets are per-cell marginal moments*. Any objective that couples months (drought duration, multi-month persistence, event-level severity) breaks the decomposition. Their §4 Conclusion explicitly identifies "explicit persistence properties... such as limiting maximum drought duration or enforcing multi-year drought conditions through constraint-based modifications" as future work, "though at the cost of increasing the computational effort to meet the targets" — this is the MOEA-FIND research question, named as future work in a same-paradigm peer paper. MOEA-FIND additionally uses Borg MM MOEA (DD-07) producing a Pareto archive that tiles the feasible drought-hazard region under the L1+anti-ideal construction (DD-11), rather than gridding a 2-D control plane.

**Code:** MIT-licensed, public on GitHub and Zenodo (10.5281/zenodo.17195525). Vendored under `external_source_code/Gozinih-Target-Synthetic-Streamflow-Generator-fea78b9/` for offline reference. Written in Python with Numba `@njit` decorators on the Kirsch core; `joblib.Parallel` over locations; SciPy `differential_evolution` per (month, location). User-facing UX is `.xlsx` inputs and a `Run.bat` driver — practitioner-flavored, not HPC-flavored.

**Detailed notes:** `manuscript/literature/notes/gozini_2026.md`.

### Herman et al. and Reed Group — MORDM Applications
The Reed group (Cornell) extensively employs synthetic streamflow generation within the Many-Objective Robust Decision Making (MORDM) framework. Synthetic generation enables: (1) preserving effects of deeply uncertain factors, (2) evaluating policies across plausible futures, (3) discovering system vulnerabilities via ensemble simulation. OpenMORDM software documents this workflow.

**Relevance:** MOEA-FIND fits naturally as the scenario generation component of MORDM, providing structured drought coverage that unstructured random ensembles lack.

---

## 1B. Scenario-Neutral / Climate Stress-Testing Lineage

This section captures the parallel methodological lineage that runs through Prudhomme, Brown, Culley, Guo, and Fowler — the climate-stress-testing tradition. These authors share substantial vocabulary with the MORDM/Reed and MOEA-FIND literature but originate from a different community (Adelaide, Melbourne, Massachusetts, Wallingford) and use overlapping but not identical terminology. Several of these papers introduce or formalize terms (exposure space, scenario-neutral, climate attribute, perturbation, failure surface, critical climate condition) that MOEA-FIND must navigate carefully. Entries are ordered chronologically; the foundational originating-citation papers (Prudhomme 2010, Brown 2012) appear first, followed by methodological extensions (Guo 2018, Guo 2017, Culley 2021) and the modern review (Fowler 2024).

### Prudhomme et al. (2010) — Scenario-Neutral Origin
**Citation:** Prudhomme, C., Wilby, R.L., Crooks, S., Kay, A.L., and Reynard, N.S. (2010). Scenario-neutral approach to climate change impact studies: Application to flood risk. *Journal of Hydrology*, 390, 198-209. doi:10.1016/j.jhydrol.2010.06.043

**Originating paper for the term "scenario-neutral."** Defines the approach (page 198 abstract): "based on sensitivity analyses of catchment responses to a plausible range of climate changes (rather than the time-varying outcome of individual scenarios), making it scenario-neutral." Frames the framework as separating two distinct elements: "the climate change projections (the hazard) from the catchment responsiveness (the vulnerability) expressed as changes in peak flows." This hazard / vulnerability separation is the conceptual bedrock of every subsequent scenario-neutral / stress-testing paper.

Originating use of the **change factor (CF) / delta change / perturbation method** (page 201): the technique of applying GCM-derived percentage or absolute changes to a reference climatology to construct perturbed time series for impact modeling. Three steps: (1) define reference climatology, (2) calculate CFs from GCM grid-box closest to target, (3) add change to reference series.

Introduces a vocabulary cluster of synonyms for safety margins: "climate change allowance," "design allowance," "headroom," and "freeboard" (page 199) — all referring to the precautionary increment built into engineering design for climate uncertainty.

Demonstrates the framework by testing the UK Government's 20% flood-flow allowance against an ensemble of 17 GCMs and 3 SRES emission scenarios at two contrasting catchments (Enrick, Roding). Constructs response surfaces in 2-D (mean precipitation change × seasonality amplitude) and overlays GCM ensemble points to assess how many projections fall outside the 20% safety margin.

**Relevance:** Foundational for vocabulary. Cite for: "scenario-neutral" (originating term), "hazard / vulnerability separation," "change factor / delta change / perturbation," "response surface" (early use), and the broader framing of bottom-up climate impact assessment as the inversion of top-down GCM-led assessment. Note that Prudhomme's response surface is in *climate-stressor* coordinates (mean precipitation × seasonality), placing this paper firmly in Fowler's "climate stressor mode" exposure space (Fowler 2024 page 2; Culley 2016 Figure 3 lineage).

### Brown et al. (2012) — Decision Scaling Origin
**Citation:** Brown, C., Ghile, Y., Laverty, M., and Li, K. (2012). Decision scaling: Linking bottom-up vulnerability analysis with climate projections in the water sector. *Water Resources Research*, 48, W09537. doi:10.1029/2011WR011212

**Originating paper for the term "decision scaling."** Definition (page 1, paragraph 5): "decision-scaling refers to the use of a decision analytic framework to reveal the scaling of climate information that is needed to best inform the decision at hand." The premise inverts the conventional process: instead of starting with GCM projections and propagating them through models to estimate impacts, decision scaling starts with decision analysis to identify climate states that favor particular decisions, then uses GCM projections only afterward to estimate the relative probabilities of those climate states.

Three-step methodology (Figure 1, pages 3-7):
1. **Identification of climate concerns, hazards, and thresholds** — through stakeholder dialogue and historical record review. Build a "decision system model" relating climate inputs to system performance.
2. **Risk discovery** — sensitivity analysis using a large stochastic input series (tens of thousands of years) to sample possible climate conditions. Develop a **climate response function** g(v_T) where Y_T = g(v_T): a surrogate model linking climate variables to performance indicators.
3. **Tailoring climate information** — use GCM ensembles to estimate probabilities for the climate states identified in step 2.

Introduces or formalizes several closely linked terms:
- **Climate state** — "the range of climate variables that favors a particular decision option" (page 3, paragraph 15). The decision-theoretic ξ_i^C in equation (1).
- **Decision threshold** — "a point where the optimal decision changes as a function of the climate conditions" (page 3, paragraph 15).
- **Climate response function** — "a surrogate model, representing the results of a series of models in a computationally efficient form that links climate variables directly to performance indicators" (page 4, paragraph 25). g(v_T) where v_T is a reduced set of influential climate statistics.
- **Climate space** — used (page 5, paragraph 21) for the multi-dimensional space of climate variables; this is the direct predecessor of Culley 2016's and Fowler 2024's "exposure space."

Critical framing claim (page 2, paragraph 5): GCM projections "represent the irreducible lower bound on the range of climate uncertainty [Stainforth et al., 2007], they should not be used to identify risks, but rather as a potential prioritization weighting on risks." This is a foundational philosophical claim of the bottom-up paradigm — GCMs constrain probability assignments, not hazard identification.

**Relevance:** Foundational for vocabulary. Cite for: "decision scaling" (originating term), "climate response function," "climate state," "decision threshold," and the inversion of the top-down workflow. Note that Brown's "climate space" is the direct predecessor of "exposure space" (Culley 2016 attribution in Fowler 2024 Table 1). The three-step decision-scaling workflow and the climate-response-function surrogate-modeling idea anchor much of the subsequent literature including the Reed-group EMODPS framing.

### Fowler et al. (2024) — Climate Stress Testing Review
**Citation:** Fowler, K.J.A., McMahon, T.A., Westra, S., Horne, A., Guillaume, J.H.A., Guo, D., Nathan, R., Maier, H.R., and John, A. (2024). Climate stress testing for water systems: Review and guide for applications. *WIREs Water*, 11(6), e1747. doi:10.1002/wat2.1747

The most current practitioner-oriented review of climate stress testing for water systems. **Most importantly for MOEA-FIND, Table 1 (pages 3-4) is an explicit glossary of 18 terms** (adaptive capacity, baseline, bottom-up, decision scaling, deep uncertainty, exposure space, failure domain, failure point, failure surface, non-climate stressor, perturbation, scenario discovery, scenario-neutral, sensitivity analysis, stochastic sequence, stress testing, stressor, top-down, vulnerability analysis) with originating-citation references for each. This is the canonical anchor for the field's vocabulary and should be the primary reference for the MOEA-FIND glossary. The paper also provides:

- A two-phase workflow proposal (minimal-scope first — perturb streamflow directly to find failure regions; expand-scope second — perturb climate with rainfall-runoff modeling and validate against phase 1) that is structurally analogous to MOEA-FIND's "search hazard space first, characterize ψ later" framing.
- A four-way taxonomy of perturbation methods: (1) simple historical scaling, (2) resampled-baseline scaling, (3) perturbing stochastic generator parameters directly, (4) manipulation of a stochastic baseline. Gozini (2026), Wheeler (2025), and MOEA-FIND all live in (3); FIND (Zaniolo 2024) lives in (4).
- A four-way taxonomy of stress-test result visualization: (i) performance vs. stressors via line/contour/colored regions; (ii) individual scenarios for in-depth discussion; (iii) summarizing failure domain (scenario discovery / PRIM / decision trees); (iv) characterizing failure boundary as a manifold (MORE/POMORE optimization, info-gap).
- Tables S1, S2 list 82 stochastic generator models classified by stationary/non-stationary, single/multi-site, and timescale (annual through daily).
- Box 2 contrasts streamflow-axes (Brown & Wilby 2012) vs. climate-axes (Mukundan et al. 2019) exposure-space designs — directly relevant to the modeling-scope decision.

**Relevance:** Top-priority citation for MOEA-FIND's vocabulary and framing. Fowler 2024 provides the modern unified language of the bottom-up/scenario-neutral/stress-testing field. The forthcoming MOEA-FIND glossary should be grounded in Fowler 2024 Table 1 and cross-reference each definition to its originating citation. Note that "exposure space" in Fowler is the canonical term; "scenario-neutral space" (Culley 2021) and "response surface" are listed as synonymous. MOEA-FIND should pick one term and stick with it; recommend "exposure space" as the dominant usage.

### Guo et al. (2018) — Inverse Approach to Rainfall Perturbation
**Citation:** Guo, D., Westra, S., and Maier, H.R. (2018). An inverse approach to perturb historical rainfall data for scenario-neutral climate impact studies. *Journal of Hydrology*, 556, 877-890. doi:10.1016/j.jhydrol.2016.03.025

(Note: Available-online date is 22 March 2016 in the journal's citation system, hence Gozini et al. 2026 cites this as "Guo et al. 2016.") Originates the **inverse approach** to scenario-neutral exposure-space generation: optimize the parameters of a parametric stochastic weather generator (WGEN) such that the resulting time series achieve user-specified target attribute values. Demonstrates evenly-distributed sampling over the exposure space at Adelaide, Australia. Compares forward approach (perturb generator parameters and accept the resulting attribute values) versus inverse approach (perturb generator parameters such that resulting attribute values match a regular target grid). This is the foundational "inverse approach" paper that Borgomeo (2015), Wheeler (2025), Gozini (2026), and MOEA-FIND all extend or reference. Six-pillar exposure-space-attribute framework: long-term mean, extremes (99th percentile), seasonality, intermittency (annual wet days, average dry-spell length), and inter-annual variability.

**Relevance:** Direct origin of the term "inverse approach" as a method category. Fowler 2024 (page 19) lists three responses to the parameter→attribute non-linearity problem: (1) accept and live with non-uniform coverage (Dubrovsky et al. 2000), (2) post-process via quantile mapping (Steinschneider & Brown 2013), or (3) **apply an inverse approach via optimization** (Guo et al. 2018). MOEA-FIND should cite Guo 2018 (alongside Borgomeo 2015) as the originating inverse-approach reference, distinct from the Borgomeo/Wheeler nonparametric-resampling lineage.

### Guo et al. (2017) — Hydro-Meteorological Attribute Sensitivity
**Citation:** Guo, D., Westra, S., and Maier, H.R. (2017). Use of a scenario-neutral approach to identify the key hydro-meteorological attributes that impact runoff from a natural catchment. *Journal of Hydrology*, 554, 317-330. doi:10.1016/j.jhydrol.2017.09.021

Extends Guo 2018's inverse approach by perturbing six hydro-meteorological attributes simultaneously (annual rainfall, winter rainfall, 99th-percentile daily rainfall, average dry-spell length, mean temperature, mean potential evapotranspiration) and applying Sobol' global sensitivity analysis to identify the dominant attributes for catchment runoff. Result: winter rainfall and annual rainfall dominate; PET-related variables matter little. Establishes that **the choice of which attributes to include in the exposure space is itself a consequential decision** — different attributes produce different sensitivity rankings.

**Relevance:** Important precedent for the argument that scenario-neutral / stress-testing analysis is sensitive to the choice of exposure-space dimensions, not just to the points sampled within. This is a transferable critique to Gozini's 2-D (mean × SD) exposure space — a higher-dimensional analysis might reveal sensitivities masked by the moment-only framing. MOEA-FIND can reuse this argument structure for the moment-vs-event-level framing.

### Culley et al. (2021) — Critical Climate Conditions Selection
**Citation:** Culley, S., Maier, H.R., Westra, S., and Bennett, B. (2021). Identifying critical climate conditions for use in scenario-neutral climate impact assessments. *Environmental Modelling and Software*, 136, 104948. doi:10.1016/j.envsoft.2020.104948

Formalizes the four-stage scenario-neutral workflow:
1. Selection of climate attributes (axes of the scenario-neutral space)
2. Development of perturbed attribute values (sampling locations in the space)
3. Generation of climate-perturbed time series
4. System performance assessment

Develops a method for selecting **critical climate attributes** — narrowing a long candidate list (15 attributes in the Lake Como case study) to a shortlist of those most relevant for the system's performance, via sparse sensitivity-analysis screening. Demonstrates that different performance objectives (flood reliability vs. irrigation deficit at Lake Como) depend on different critical attributes, so the axes-selection stage must be tailored to the objective. Uses the term "scenario-neutral space" interchangeably with "exposure space" and "response surface."

**Relevance:** Critical for the manuscript's framing of stage (i) — selection of axes for the exposure / hazard space. MOEA-FIND's use of three SSI-derived event-level characteristics (mean event severity, mean cumulative deficit, time-in-drought fraction) as the K=3 production hazard space (DD-04) is a stage-(i) decision in Culley 2021 terms. The Culley 2021 framework gives MOEA-FIND a credible methodological story for *why* these three axes were chosen — they are the system-relevant emergent characteristics for drought-driven failure modes — and a precedent for arguing that the axes-selection decision is itself a decision-relevant methodological fork.

**Synonym alert:** Culley 2021 uses "scenario-neutral space" as a *new* synonym for what Brown 2012 called "climate space," what Quinn 2018 calls "exposure space," and what some literature calls "response surface." Fowler 2024 Table 1 lists "exposure space" as the canonical term. The MOEA-FIND glossary must catalogue these synonyms and explain the chosen usage.

### Cross-paper terminology observations (Fowler, Guo, Culley)

Reading the four papers as a set, several inter-author synonym/conflict patterns emerge that the glossary will need to address:

1. **Exposure space / scenario-neutral space / response surface / climate space.** Fowler 2024 standardizes on "exposure space" with originating reference to Brown 2012 ("climate space"). Culley 2021 introduces "scenario-neutral space." All four terms refer to the same multi-dimensional space of stressors over which performance is mapped. **Recommended MOEA-FIND choice:** "exposure space," but acknowledge Culley's "scenario-neutral space" as field-current.

2. **Climate attribute / climate variable / hydrologic stressor / hydro-meteorological attribute / stressor.** Guo 2017 says "hydro-meteorological attribute" (a specific statistic of a variable, e.g., 99th percentile of daily rainfall). Culley 2021 says "climate attribute" for the same concept. Fowler 2024 reserves "stressor" for the axis of the exposure space and "variable" for the underlying physical quantity. **Recommended MOEA-FIND choice:** "stressor" for an axis of the exposure space; "drought characteristic" for a hazard-space axis.

3. **Bottom-up vs. scenario-neutral vs. stress testing.** Fowler 2024 §2 carefully distinguishes: bottom-up is the broader paradigm; stress testing is a specific type of bottom-up sensitivity analysis with focus on identifying failure modes; scenario-neutral is a specific name for stress-testing methods that operate independently of GCM scenarios (Prudhomme 2010). **Recommended MOEA-FIND choice:** "stress testing" for the workflow; "bottom-up" for the paradigm; "scenario-neutral" only when contrasting against top-down GCM-driven approaches.

4. **Failure point / failure domain / failure surface.** Fowler 2024 codifies these (after Bryant & Lempert 2010, Guillaume et al. 2016): failure point is one combination of stressors that yields failure; failure domain is the contiguous region of the exposure space in which failure occurs; failure surface is the boundary between failure and non-failure regions. The Culley 2021 phrase "regions of success and failure" maps to "non-failure / failure domains."

5. **Inverse approach.** Originates with Guo 2018 in the parametric weather-generator context. Borgomeo 2015 used the same method with a non-parametric generator without using the term "inverse." Fowler 2024 (page 19) lists "inverse approach" as one of three responses to the parameter→attribute non-linearity. Gozini 2026 adopts the term explicitly. **Recommended MOEA-FIND choice:** "inverse approach" is good shared vocabulary; cite Guo 2018 alongside Borgomeo 2015.

---

## 2. Kirsch-Nowak Modified Fractional Gaussian Noise Methods

### Kirsch et al. (2013)
**Citation:** Kirsch, B.R., Characklis, G.W., and Zeff, H.B. (2013). Evaluating the impact of alternative hydro-climate scenarios on transfer agreements: Practical improvement for generating synthetic streamflows. Journal of Water Resources Planning and Management, 139(4), 396-406.

Modified fractional Gaussian noise with uncorrelated monthly bootstrap resampling followed by Cholesky correlation restoration. Operates via monthly resampling with independent draws from historical records, then imposing historical monthly correlation structure through Cholesky decomposition of the twelve-by-twelve monthly correlation matrix in Gaussian space. Addresses Fractional Gaussian Noise limitations while preserving seasonal correlation. Widely implemented.

**Relevance:** Core generation mechanism for MOEA-FIND's nonparametric pathway.

### Nowak et al. (2010)
**Citation:** Nowak, K., et al. (2010). A nonparametric stochastic approach for multisite disaggregation of annual to daily streamflow. Water Resources Research, 46, W08529.

Combines nonhomogeneous Markov chain with K-nearest neighbor (K-NN) time series bootstrap conditioned on hydrologic state. Reproduces daily data at multiple sites while preserving summability, continuity, and cross-site dependencies.

**Relevance:** K-NN conditioning could inform how MOEA-FIND selects bootstrap blocks.

### Block Size Selection
Block size critically affects bootstrap validity. Methods in the literature:
- **Hall, Horowitz, and Jing (HHJ):** Cross-validation to minimize MSE across block lengths
- **Nonparametric Plug-In (NPPI):** Estimates theoretically optimal block length via resampling
- **Variable-Length Block (VLB):** Weighted fragment methods for disaggregation with perturbation; generates flows outside historical record

**Relevance:** MOEA-FIND's block size choice (DD-03) should be informed by these methods. VLB is particularly interesting as it can extend beyond historical range.

---

## 3. Parametric Streamflow Generation with Copulas

### Kappa Distribution (Svensson et al. 2017)
**Citation:** Svensson, C., et al. (2017). Statistical distributions for monthly aggregations of precipitation and streamflow in drought indicator applications. Water Resources Research, 53(12).

Four-parameter kappa distribution provides excellent fit for monthly streamflow across physiographic regions. Flexible tail behavior, bounded below at zero. Superior for low-flow volumes critical to drought characterization.

**Relevance:** Recommended marginal distribution for MOEA-FIND's parametric CDF pathway.

### Vine Copulas for Streamflow
Multiple recent papers establish vine copulas as state-of-the-art for temporal/spatial dependence:

- **PAR(p)-vine copula models** (Brechmann et al., 2017, Stoch. Environ. Res. Risk Assess.): Combine periodic autoregression with vine copulas for stochastic streamflow scenario generation.
- **D-vine copulas with lag-2** (various, 2019-2023): Capture temporal autocorrelation at lag-1 and lag-2 simultaneously. Mixed D-vine conditional quantile model for monthly simulation (Li et al. 2023, J. Hydro-environ. Res.).
- **Periodic spatial vine copulas** (2019, Water Resour. Manage.): Account for seasonally varying dependence structure.
- **Spatio-temporal vine copulas** (2018, Water): Synthesize continuous precipitation/streamflow with vine copula dependence.

**Relevance:** When MOEA-FIND uses parametric marginals (kappa CDFs), vine copulas impose the temporal and spatial correlation structure.

### Papalexiou and Koutsoyiannis — Stochastic Hydrology Advances
Foundational contributions: (1) long-range dependence and second-order structure identification, (2) frequency-domain simulation via phase randomization, (3) CoSMoS framework for flexible distribution modeling with copula-based intermittency (Papalexiou, 2022, WRR). Phase randomization + kappa distribution allows extrapolation to unobserved extremes while preserving spectral properties.

**Relevance:** Phase randomization is a potential third pathway for MOEA-FIND (beyond bootstrap and parametric CDF). The CoSMoS framework may offer a unified generation approach.

### KDE-Based Approaches
Mixed copula-KDE framework (2025, Stoch. Environ. Res. Risk Assess.) demonstrates nonparametric smoothing + copula dependence. Extends support beyond observed range via kernel bandwidth without parametric assumptions.

**Relevance:** A middle-ground option for MOEA-FIND between pure bootstrap and parametric CDF.

---

## 4. Multi-Objective Optimization for Scenario Design

### Borg MOEA
**Citation:** Hadka, D. and Reed, P. (2013). Borg: An auto-adaptive many-objective evolutionary computing framework. Evolutionary Computation, 21(2).

Three key mechanisms make Borg ideal for MOEA-FIND:
1. **Epsilon-dominance archiving:** Maintains diverse approximate solutions while controlling archive size
2. **Epsilon-progress restart triggering:** Escapes premature convergence
3. **Adaptive operator selection:** Six variational operators (SBX, DE, PCX, UNDX, SPX, PM, UM) adapt to problem-specific landscape

Demonstrated on 33 benchmark instances (DTLZ, WFG, CEC 2009). Particular success in environmental applications.

**Relevance:** Epsilon-dominance archiving is the mechanism that produces near-uniform coverage on the Manhattan-norm hyperplane. No other MOEA has this specific combination.

### Manhattan/L1 Norm in Multi-Objective Optimization
The "minimum Manhattan distance" (MMD) approach appears in solution selection literature: select solutions corresponding to minimal L1 distance from an ideal vector. In MOEA-FIND, the L1 norm is used as an *objective* rather than a selection criterion, forcing all Pareto-optimal solutions onto a hyperplane where J_1 + J_2 + ... + J_k + J_{k+1} = constant. This is a novel application of the L1 norm in evolutionary computation.

### Reference-Point MOEAs (NSGA-III, MOEA/D)
- **NSGA-III** (Deb and Jain, 2014, IEEE Trans. Evol. Comput.): Uses pre-defined reference points for uniform Pareto front coverage
- **MOEA/D** (Zhang and Li, 2007, IEEE Trans. Evol. Comput.): Decomposes into scalar subproblems

These achieve coverage through different mechanisms than epsilon-dominance. Not recommended for MOEA-FIND because they lack the specific hyperplane-tiling property that the Manhattan norm + epsilon-dominance provides.

### MOEA for Scenario Generation (Novel Application)
No published work uses MOEAs specifically for structured scenario generation where drought characteristics are objectives. This is the core novelty of MOEA-FIND. The closest precedent is greedy ensemble selection from Pareto fronts in the Reed group's work, but applied to policy selection, not scenario generation.

---

## 5. Drought Characterization and Metrics

### Standard Indices
- **SPI** (McKee et al.): Precipitation only. Underestimates drought frequency by 18-27% in drylands.
- **SPEI** (Vicente-Serrano et al.): Includes temperature via reference ET. Better for warming scenarios.
- **SSI** (Standardized Streamflow Index): Applied directly to streamflows. Captures hydrological drought. Used by Zaniolo (FIND) for drought control.

### Severity-Duration-Frequency (SDF) Curves
Bivariate frequency analysis (often copula-based) produces SDF curves showing joint distributions of drought magnitude and duration. Standard tools for reservoir design and drought mitigation. DSDF curves inform the range of drought characteristics that MOEA-FIND should target.

### Threshold Selection
- **Fixed percentile:** Q80 (80th percentile of flow duration curve) or P20 (20th percentile of monthly flow). Year-round constant.
- **Variable daily threshold (VDT):** 95th percentile computed via 31-day moving window. Accounts for seasonal variation.
- **System-specific:** E.g., FFMP drought trigger levels for NYC DRB.
- SSI convention: drought onset at SSI less than or equal to -0.84 (20th percentile).

**Relevance:** MOEA-FIND's drought metrics must use consistent, defensible thresholds. Recommend P20 monthly flow as default, with system-specific thresholds for case studies.

### Multi-Dimensional Characterization
Drought is not a single phenomenon. Comprehensive monitoring integrates meteorological (weeks), agricultural (months), and hydrological (months-years) timescales. MOEA-FIND focuses on hydrological drought (streamflow-based), with event-level metrics: duration, peak intensity, cumulative severity, and (for longer traces) frequency.

---

## 6. Policy Re-Evaluation, Robustness, and Scenario Discovery

### EMODPS Framework
The Reed group's EMODPS pairs multi-objective evolutionary optimization of operating policies with visual analytics for tradeoff discovery. Generated scenarios stress-test optimized policies. This is the primary consumer of MOEA-FIND output.

### Bryant and Lempert (2010) — Scenario Discovery Origin
**Citation:** Bryant, B.P. and Lempert, R.J. (2010). Thinking inside the box: A participatory, computer-assisted approach to scenario discovery. *Technological Forecasting and Social Change*, 77, 34-49. doi:10.1016/j.techfore.2009.08.002

**Originating paper for the term "scenario discovery."** Definition (page 1, abstract): "scenario discovery defines scenarios as a set of plausible future states of the world that represent vulnerabilities of proposed policies, that is, cases where a policy fails to meet its performance goals." The approach uses statistical or data-mining algorithms applied to large databases of simulation-model results to identify "easy-to-interpret combinations of uncertain model input parameters that are highly predictive of these policy-relevant cases."

Formal four-step workflow (Figure 1, page 4):
1. Generate data from a simulation model (typically by Latin Hypercube sampling over uncertain inputs).
2. Identify candidate scenarios using PRIM (Patient Rule Induction Method, originating from Friedman & Fisher 1999) or CART (Classification and Regression Trees).
3. Assess scenarios with diagnostics (resampling test, quasi-p-value test).
4. Choose scenarios.

Introduces or formalizes critical vocabulary:
- **Vulnerability** (page 2, paragraph following Fig. 1 cite): "states of the world where a proposed policy may fail to meet its performance goals as well as those where a policy's performance deviates significantly from the optimum policy in that state of the world." Footnote 4: "The former represents an absolute performance measure and the latter a regret measure."
- **States of the world** — used in the formal sense of Lempert decision theory: points in M-dimensional space of uncertain model input parameters.
- **Cases of interest** (I_s) — set of cases where some policy-relevant criterion crosses a threshold; the targets of scenario discovery. Fowler 2024 Table 1 cites this exact term as the origin for "failure point."
- **Box / box set** — a multi-dimensional region defined by constraints {a_j ≤ x_j ≤ b_j, j ∈ L_k} on a subset of input parameters; each box is interpreted as a scenario, the unconstrained parameters become the "key driving forces."
- **Coverage / density / interpretability** — three measures of merit for scenario quality:
    - *Coverage* = fraction of policy-relevant cases captured by the scenarios (analogous to recall/sensitivity in classification).
    - *Density* = fraction of cases inside the box that are policy-relevant (analogous to precision/positive predictive value).
    - *Interpretability* = qualitative measure approximated by counting parameters constrained per box; "highly interpretable box set should consist of on the order of three or four boxes, each with on the order of two or three constrained parameters" (page 5).
- **Key driving forces** — input parameters that the algorithm constrains in defining a box; explicitly analogous to the "key driving forces" of the intuitive logics scenario school (Schwartz 1991).
- **Deep (Knightian) uncertainty** (footnote 1, page 2): "the condition where parties to the decision do not know or do not agree on the system model relating actions to consequences or the prior probability distributions for input parameters to these system models."

PRIM operation summary: a bump-hunting algorithm that iteratively peels off the worst-performing slice of the input space (steps 1–N "peeling") then refines by re-pasting (steps N+1 onwards). Generates a coverage-density tradeoff curve. Paper provides software (sdtoolkit R package) and is the canonical implementation reference for PRIM in scenario discovery.

**Relevance:** Foundational citation for scenario discovery. Cite for: "scenario discovery" (originating term), "vulnerability" (formal definition), "cases of interest," "coverage," "density," "interpretability" (originating definitions), the formal four-step workflow, and the canonical PRIM implementation. The Hadjimichael et al. (2020 onward) and Bonham (2024, 2025) literature builds directly on this foundation. MOEA-FIND's coverage / density / interpretability framing in the Pareto archive sense is a deliberate reuse of this vocabulary at a different stage of the workflow (scenario generation rather than scenario discovery) — the manuscript should acknowledge the lineage to avoid reader confusion.

**Synonym alert:** "Cases of interest" (Bryant & Lempert) ≈ "failure points" (Fowler 2024 Table 1). "States of the world" (Bryant & Lempert) ≈ "states of the world" (Hadjimichael 2020) but in different sampling contexts (parameter LHS vs. SOW + inner stochastic). "Box" / "box set" terminology is specific to Bryant & Lempert / PRIM and is not used in the Fowler 2024 / Reed-group lineage; that literature uses "scenario" directly.

### BART for Scenario Discovery
**Citation:** Chipman, H.A., et al. (2010). BART: Bayesian additive regression trees. Annals of Applied Statistics, 4(1).
**Application:** Water Programming blog (2024) documents BART as the preferred modern method for scenario discovery, replacing PRIM. Offers probabilistic uncertainty quantification and superior scalability.

**Relevance:** MOEA-FIND ensembles feed into BART-based vulnerability analysis. Structured coverage ensures BART has diverse training data across drought characteristic space.

### Hadjimichael et al. — Stakeholder-Contingent Robustness
**Citation:** Hadjimichael, A., et al. (2020). Defining Robustness, Vulnerabilities, and Consequential Scenarios for Diverse Stakeholder Interests in Institutionally Complex River Basins. Earth's Future, 8(6).

Robustness is stakeholder-specific: the same drought affects different stakeholders differently. Satisficing thresholds (minimum acceptable performance) vary by stakeholder. MOEA-FIND's uniform drought space coverage ensures all stakeholders encounter both success and failure conditions.

### Sunkara et al. (2023)
**Citation:** Sunkara, et al. (2023). How Should Diverse Stakeholder Interests Shape Evaluations of Complex Water Resources Systems Robustness When Confronting Deeply Uncertain Changes? Earth's Future, 11(7).

Extends Hadjimichael's framework. BART scenario discovery identifies factors driving consequential scenarios that differ by robustness definition. MOEA-FIND ensembles are ideally suited because structured coverage ensures unbiased discovery.

### Bonham et al. (2024) — Subsampling and Space-Filling Metrics
**Citation:** Bonham, N., Kasprzyk, J., Zagona, E., and Rajagopalan, B. (2024). Subsampling and space-filling metrics to test ensemble size for robustness analysis with a demonstration in the Colorado River Basin. Environ. Model. Softw., 172, 105933.

Closest published precedent to the library-and-subsample baseline. Uses conditioned LHS (cLHS) to subsample from a pre-generated 500-scenario library (streamflow from historical, paleo, CMIP-3; demand and initial storage from LHS). Evaluates subset quality via space-filling metrics: MST mean edge length, MST standard deviation, and minimum pairwise distance. MSTmean explains 77-91% of variance in robustness ranking accuracy (Kendall tau-b). Satisficing metrics require 50-300 scenarios; regret metrics require 400+.

**Key distinction from MOEA-FIND:** Bonham et al. subsample in uncertainty *input* space (streamflow features, demand, storage levels) for computational efficiency (fewer simulations to preserve rankings). MOEA-FIND targets structured coverage of drought *outcome* space (duration, severity, frequency) for vulnerability discovery. Input-space coverage does not guarantee outcome-space coverage, especially when the mapping from inputs to drought characteristics is nonlinear.

---

## 7. Coverage Metrics and Sampling Theory

### Latin Hypercube Sampling (LHS)
Standard space-filling design. Partitions parameter space into n^d subspaces, placing exactly one point per univariate marginal interval. Reduces clustering relative to random sampling. Widely used in environmental modeling uncertainty analysis.

### Sobol Sequences
Deterministic quasi-random low-discrepancy sequences. Faster convergence than LHS in sensitivity analyses (PLOS Computational Biology comparison). Inherently extensible. Superior space-filling properties in moderate dimensions.

### Discrepancy Measures
- **Star discrepancy (L_infinity):** Maximum deviation from uniform distribution. Theoretical bounds via Koksma-Hlawka inequality.
- **L2 discrepancy:** Mean squared deviation. More practical for empirical comparison.
- **Morris SU (Sampling Uniformity):** Specific to sensitivity analysis sampling designs.

**Relevance:** MOEA-FIND's coverage quality should be evaluated using L2 discrepancy and compared to LHS and Sobol baselines. This provides the quantitative backbone for RQ1 and RQ4.

### Pareto Front as Space-Filling Design
**Citation:** IEEE CEC (2019) work on sampling reference points on Pareto fronts of benchmark problems.

Uniformly distributing K solutions across irregularly shaped Pareto fronts is an active research problem. The SPREAD algorithm uses diffusion models for iterative refinement. MOEA-FIND's Manhattan norm trick sidesteps this by forcing the front onto a regular hyperplane where epsilon-dominance naturally provides uniformity.

**Relevance:** This is the theoretical justification for why MOEA-FIND produces structured coverage. The Pareto front IS a space-filling design in drought characteristic space, with coverage quality controlled by epsilon values.

---

## Summary of Key References by Experimental Phase

### POC (Experiments 1.1-1.3)
- Hadka and Reed (2013) — Borg MOEA
- Deb and Jain (2014) — NSGA-III comparison
- LHS/Sobol theory for coverage baselines

### Single-Site Bootstrap (Experiments 2.1-2.4)
- Kirsch et al. (2013) — Block bootstrap
- Borgomeo et al. (2015) — Search-based generation
- Zaniolo et al. (2024) — FIND comparison
- Svensson et al. (2017) — Drought threshold/metrics

### Parametric CDF (Experiments 3.1-3.3)
- Svensson et al. (2017) — Kappa distribution
- PAR(p)-vine copula (Brechmann et al. 2017)
- Papalexiou and Koutsoyiannis — Phase randomization
- Mixed copula-KDE (2025)

### Multi-Site DRB (Experiments 4.1-4.3)
- Wheeler et al. (2025) — Multi-site generation
- Nowak et al. (2010) — Multi-site disaggregation
- Hadjimichael et al. (2020) — Stakeholder robustness
- BART (Chipman et al. 2010; Water Programming 2024)

---

## 8. Extraction round 2026-04-14: fifteen additional papers

Added 2026-04-14 after two parallel opus extraction agents read every
PDF in `manuscript/literature/` that had not yet been extracted.
Quotations are verbatim with page numbers. Source: Phase 6.2 and 6.3
progress log entries in `C:\Users\tjame\.claude\plans\staged-painting-dongarra.md`.

### 8.1 Scenario discovery and robustness tradition

**Hadjimichael, Quinn, Wilson, Reed, Basdekas, Yates, Garrison (2020).**
Defining robustness, vulnerabilities, and consequential scenarios for
diverse stakeholder interests in institutionally complex river basins.
*Earth's Future* 7, e2020EF001503.

- One-sentence contribution (p. 1): "advances previous robustness
  evaluation methods by formalizing an a posteriori exploration of
  alternative definitions of robustness tailored to each user's
  unique context".
- Quotable definition (p. 1): "In bottom-up methods, the identification
  of consequential scenarios (also known as scenario discovery)
  presupposes a metric-based classification of the acceptability of
  system performance, most commonly by exploiting one or more candidate
  measures of robustness."
- Satisficing metric definition (p. 11): "The domain criterion
  satisficing metric (Starr, 1962) has been widely used in the
  robustness literature and quantifies the fraction of scenarios in
  which a desired performance is met."
- Sampling design (p. 6): "A Latin Hypercube sample (McKay et al., 1979)
  of 1,000 parameter combinations (SOWs) is generated across the ranges
  given in Table 1 assuming uniformity and independence."
- Terms defined: consequential scenarios (p. 1), satisficing threshold
  (p. 11), states of the world / SOWs (p. 6), domain criterion
  satisficing metric (p. 11), Patient Rule Induction Method (p. 11),
  a posteriori robustness (p. 1).
- Sampling claim: operates exclusively in parameter (input) space via
  LHS over 14 deeply uncertain factors. No claim of drought hazard
  space uniformity.
- Figure template: Figure 3 (SOW-to-performance pipeline diagram) and
  Figure 7 (PRIM factor maps with robustness thresholds) are direct
  templates for MOEA-FIND Figures 1 and 7 respectively.

**Hadjimichael, Quinn, Reed (2020).** Advancing diagnostic model
evaluation to better understand water shortage mechanisms in
institutionally complex river basins. *Water Resources Research* 56,
e2020WR028079.

- Contribution (p. 1): "a diagnostic evaluation framework that pairs
  exploratory modeling with global sensitivity analysis".
- Behaviour-space quotation (p. 2): "To be effective in producing a rich
  enough picture of a complex model's behavior space... exploratory
  modeling must examine a very large and diverse suite of model
  simulation runs that capture important interactions and mechanisms
  leading to consequences of interest."
- Sampling (PDF p. 7): 1,000 LHS SOWs in a 14-dim HMM hyperparameter
  plus human-system multiplier parameter space, 10 HMM realizations
  per SOW.
- Terms: behavior space (p. 2), consequence-oriented sensitivity
  analysis (p. 2), magnitude-varying sensitivity analysis (p. 1).

**Hadjimichael, Reed, Quinn (2020).** Navigating deeply uncertain
tradeoffs in harvested predator-prey systems. *Complexity* 2020,
Article ID 4170453.

- Robustness framing (p. 8): "decision makers may consider identifying
  potential policies that continue to perform satisfactorily when
  operated under a broad range of alternative system characteristics".
- Domain criterion satisficing measure (p. 8): "quantifies the fraction
  of potential SOWs, in which a solution meets a desired performance".
- Sampling (p. 8): 4,000 LHS SOWs in a 9-dimensional deeply uncertain
  parameter space, with infeasible SOWs (species extinction) pruned
  ex post. "Such SOWs were omitted from the robustness analysis... so
  as to only evaluate them in contexts where the choice of strategy
  actually matters".
- Figure templates: Figure 6 (3D scatter of parametric space coloured
  by collapse) and Figure 7 (parallel-coordinate tradeoffs).

**Hadjimichael, Reed, Quinn, Vernon, Thurber (2024).** Scenario
storyline discovery for planning in multi-actor human-natural systems
confronting change. *Earth's Future* 12, e2023EF004252.

- FRNSIC framing (p. 1): "scenario discovery framework that addresses
  these challenges by organizing and investigating consequential
  scenarios using hierarchical classification of diverse outcomes
  across actors, sectors, and scales".
- Scenario discovery definition (p. 2): "these approaches focus on
  the exploration of large ensembles of possible futures and the a
  posteriori identification of consequential scenarios".
- Key motivational quote for MOEA-FIND (p. 2): "Consequential scenarios
  identified in this manner are entirely dependent on choices analysts
  make with regard to how they aggregate across SOWs." Cited in
  MOEA-FIND §4.5 as the central DMDU limitation that hazard-space
  sampling addresses.
- Storyline definition (p. 2, citing Shepherd 2018): "a physically
  self-consistent unfolding of plausible future events".

Preprint companion (Hadjimichael et al. 2023, ESSOAr): same FRNSIC
framework, with the quotable "large ensemble exploratory modeling...
creates a tension: on one hand... we need to create large ensembles...
On the other hand, each additional dimension... makes the results...
more difficult to convey actionable insights" (PDF p. 4).

**Moallemi, Zare, Reed, Elsawah, Ryan, Bryan (2020).** Structuring
and evaluating decision support processes to enhance the robustness
of complex human-natural systems. *Environmental Modelling and
Software* 123, 104551.

- THE central vocabulary anchor for MOEA-FIND (p. 1): "Exploratory
  modelling enables robust decision-making by exploring the
  implications of varying decision assumptions (within a decision
  space) and scenarios (within an uncertainty space) in terms of
  performance (or robustness) measures (within an outcome space),
  using a series of computational experiments." The three-space
  taxonomy (decision / uncertainty / outcome) is the frame in which
  MOEA-FIND's hazard-space coverage claim is legible.
- A posteriori framing (p. 2): "candidate decisions, key uncertainties,
  their conflicting objectives, and assumed preferences are discovered
  from a diverse search... This approach is mainly driven by the mantra
  'generate then choose within the broader information context.'"
- Terms: decision space / uncertainty space / outcome space (p. 1), a
  posteriori analysis (p. 2), decision fork (p. 1).
- MOEA-FIND positioning: the three-space taxonomy is cited in §1.3
  and §4.5 of the manuscript.

**Lau, Reed, Gold (2023).** Evaluating implementation uncertainties
and defining safe operating spaces for deeply uncertain cooperative
water resource systems. *Water Resources Research* 59, e2023WR034841.

- Implementation uncertainty (p. 2): "deviations in how utilities
  operationalize their collective and individual action policy
  pathways' rule systems... how much variation in decision variables
  can be tolerated for policy actions to retain acceptable
  performance".
- Safe operating space (p. 2): "identifying... tolerable windows of
  decision deviations where they can individually and collectively
  maintain acceptable levels of performance".
- Sampling: decision-variable space, not drought hazard space.
- Figure template: SOS boundary plots and parallel-axis deviation
  plots, transferable to MOEA-FIND SI figures.

### 8.2 Synthetic generation and MOEA foundations

**Herman, Zeff, Lamontagne, Reed, Characklis (2016).** Synthetic
drought scenario generation to support bottom-up water supply
vulnerability assessments. *Journal of Water Resources Planning and
Management* 142(11), 04016050.

- Three-challenge framework for drought ensemble design (p. 04016050-1):
  "the framework for generating synthetic streamflow scenarios must be
  sufficiently flexible to create extreme drought conditions not
  observed in the historical record, to support long-term planning.
  Second, it must approximately maintain historical statistical
  properties of streamflow... Finally, to remain consistent with the
  bottom-up approach, it must provide relatively simple parameters to
  adjust the severity of scenarios in an exploratory manner".
- Method basis (p. 04016050-2): "this study adapts the streamflow
  generation method developed by Kirsch et al. (2013) using a weighted
  bootstrap resampling".
- Bottom-up positioning (p. 04016050-1): "Bottom-up methods simulate a
  larger sampling of plausible scenarios to identify those that produce
  consequential effects on the system".
- Figure template: Figure 3c severity-frequency curve with named
  historical droughts, directly adoptable for the MOEA-FIND Figure 5
  hydrology demonstration.

**Quinn, Reed, Giuliani, Castelletti, Oyler, Nicholas (2018).**
Exploring how changing monsoonal dynamics and human pressures
challenge multireservoir management. *Water Resources Research* 54,
4638 to 4662.

- Bottom-up exploratory modelling (p. 4638): "Dessai et al. (2009)
  highlight the deficiencies of classical 'predict-then-act' risk-based
  assessments... Instead, they advocate for 'bottom-up' approaches in
  which exploratory modeling techniques (Bankes, 1993) are employed
  to discover robust strategies".
- Kirsch-Nowak generator citation chain (p. 4641): "Monthly flows on
  each tributary are generated synthetically using Cholesky
  decomposition of resampled historical monthly flows and then
  disaggregated to daily flows using a nearest neighbor approach
  introduced by Nowak et al. (2010)... The Cholesky decomposition
  preserves temporal correlation within each flow series, while
  resampling the same historical flows at each site preserves spatial
  correlation across them (Herman et al., 2016; Kirsch et al., 2013)."
- Adequacy-of-coverage quote (p. 4639): "it is important to ensure
  that the scenarios sampled in the exploratory analysis adequately
  capture the system dynamics that might emerge under alternative
  climatic and socioeconomic futures."

**Quinn, Reed, Giuliani, Castelletti (2019).** What is controlling our
control rules? *Water Resources Research* 55, 5962 to 5984. The paper
does not address ensemble construction critique. Its relevance to
MOEA-FIND is the convention of simulating over 100,000 years of
synthetic streamflows for stable performance estimates and the EMODPS
framing (p. 5963).

**Quinn, Hadjimichael, Reed, Steinschneider (2020).** Can exploratory
modeling of water scarcity vulnerabilities and robustness be
scenario neutral? *Earth's Future* 8, e2020EF001650.

- Scenario-neutral definition (p. 2 of 25): "scenarios are generated
  by sampling possible values of uncertain climatic and socioeconomic
  variables independently and uniformly over expansive ranges".
- Saltelli citation on neutrality (p. 3 of 25): "As noted by Saltelli
  et al. (2020), 'the technique is never neutral'; rather 'the choice
  of the methodology conditions the narrative produced by an
  analysis.'"
- Direct MOEA-FIND motivation quote (p. 3 of 25): "we argue that
  performing such sensitivity analyses may be equally important in
  the vulnerability analysis, so planners should assess the sensitivity
  and robustness of alternative water management policies using
  multiple scenario designs and assumed probability distributions."
- Figure 1 (p. 4 of 25): success/failure region overlay with
  multi-source sample clouds is a direct template for a MOEA-FIND
  figure comparing ensemble provenance.

**Bonham, Kasprzyk, Zagona (2025).** Taxonomy of purposes, methods,
and recommendations for vulnerability analysis. *Environmental
Modelling and Software* 183, 106269.

- Five purposes of vulnerability analysis (p. 2): "scoping, policy
  creation, policy comparison, negotiation and compromise, and
  monitoring and adaptation". MOEA-FIND ensembles support the first
  two purposes directly.
- Space-filling design definition (p. 3, §2.1.2): "A space-filling
  design of the SOW ensemble uses statistical designs of experiments
  to comprehensively model the relationship between SOWs and system
  performance... an algorithm creates a SOW ensemble that maximizes
  coverage of the uncertainty space".
- Why space-filling matters (p. 3): "Space-filling designs... provide
  a uniform, continuous sampling across each uncertain factor, which
  enables the factor mapping step to determine the boundary between
  different decision-relevant outcomes more accurately."
- Ensemble bias quote, MOEA-FIND motivation (§2.1.2 p. 3): "SOW
  ensembles derived from climate projections tend to contain
  relatively more 'moderate' SOWs and relatively fewer 'challenging'
  SOWs (Reis and Shortridge, 2020)."
- Hybrid subsampling precedent (§2.1.2): describes Bonham et al.
  (2022) cLHS subsampling of 500 SOWs from a ~2M SOW ensemble.
- Figure template: Figure 5 (PRIM versus tree-based factor maps) is
  the direct template for MOEA-FIND Figure 7. Figure 8 (sampling-gap
  plot for full-factorial vs space-filling) is a near-direct
  template for MOEA-FIND Figure 1.

**Herman, Quinn, Steinschneider, Giuliani, Fletcher (2020).** Climate
adaptation as a control problem. *Water Resources Research* 56,
e2019WR025502.

- Research Gap quote for MOEA-FIND Introduction (p. 13 of 32):
  "Process-based insight for synthetic generation: Many studies
  characterize uncertainty in future scenarios based on coarse-timescale
  GCM statistics, such as annual precipitation. There are opportunities
  to leverage insight into climate processes and model errors to
  inform finer-scale uncertainty characterization in synthetic
  scenarios."
- Quotable expansion (p. 14 of 32): "Leveraging such process-based
  insights, there is an opportunity for new synthetic scenario
  generation methods to create plausible dynamic projections... In
  principle, any parameter of a stochastic generator can be perturbed
  based on climate information, for example, streamflow seasonality...
  or the frequency and severity of drought events (Herman et al.,
  2016)."
- Terms: dynamic planning (p. 1), exogenous vs endogenous uncertainty
  (p. 1), cascade of uncertainty (p. 1, via Wilby and Dessai 2010).

**Hadka, Reed (2015).** Large-scale parallelization of the Borg
multiobjective evolutionary algorithm to enhance the management of
complex environmental systems. *Environmental Modelling and
Software* 69, 353 to 369.

- Seed convention (p. 361): "Each run was repeated 50 times with
  different initial random seeds so that the expected search quality
  and its deviation can be calculated." This is the Reed-group
  convention; MOEA-FIND uses the multi-seed replication convention
  without committing to exactly 50 seeds.
- Epsilon vector reporting format (Table 2, p. 361): per-objective
  precisions are tabulated alongside objective names and direction.
  Adopted as the MOEA-FIND §2.7 epsilon table format.
- Runtime diagnostics (p. 362): "runtime data is collected every
  10,000 NFE and stored in the database".

**Reed, Hadka, Herman, Kasprzyk, Kollat (2013).** Evolutionary
multiobjective optimization in water resources: the past, present,
and future. *Advances in Water Resources* 51, 438 to 456.

- Epsilon-dominance citation anchor (p. 441): "Laumanns et al. [25]
  shows that epsilon-dominance archiving satisfies elite preservation
  and guarantees a bounded archive size that varies with the decision
  maker's precision goals for each objective. The epsilon-dominance
  concept replaces point-based non-domination sorting with a
  grid-based sort."
- Performance metric triad (p. 443): "our results and analysis will
  focus on three performance measures: (1) generational distance
  [76], (2) additive epsilon-indicator [74], and (3) hypervolume
  [74]." MOEA-FIND SI reports hypervolume alongside the new
  coverage metrics.
- LHS of parameterisations (p. 443): "Latin hypercube sampling (LHS)
  was used to sample the full feasible parameter space for each
  algorithm... the runs of the LHS parameter block were replicated
  for 50 random number generator seeds."
- Control map convention (p. 443): population size versus NFE
  projection plots of performance; a template for an optional
  MOEA-FIND SI epsilon and NFE sensitivity figure.

### 8.3 Overall observations for MOEA-FIND positioning

1. Every paper in this extraction samples in parameter or
   decision-variable space. None claims uniformity of the drought
   hazard space. MOEA-FIND's hazard-space sampling claim is factually
   defensible against this corpus.
2. Moallemi 2020's three-space taxonomy (decision / uncertainty /
   outcome) is the single most useful vocabulary anchor and is cited
   in §1.3 and §4.5 of the manuscript.
3. Herman 2020's Research Gap 1 (p. 13 of 32) is the strongest single
   motivation citation and is cited in §1.1.
4. Hadjimichael 2024 p. 2 (consequential scenarios depend on
   aggregation choices) is the strongest citation for the argument
   that structuring the ensemble at generation time resolves an
   aggregation-arbitrariness problem. Cited in §4.5.
5. Hadjimichael 2020a Figure 7 (PRIM factor maps) and Bonham 2025
   Figure 5a (PRIM vs trees) are the direct templates for MOEA-FIND
   Figure 7. PRIM is the canonical scenario-discovery visual in this
   lineage; the earlier plan's gradient-boosted-tree decision boundary
   was the wrong genre and has been replaced.
6. Hadjimichael 2020a Figure 3 (SOW-to-performance pipeline) is the
   direct template for MOEA-FIND Figure 1 (parameter-space to
   hazard-space contrast).
7. Kirsch-Nowak citation chain: Kirsch 2013 + Nowak 2010 + Herman 2016
   is the Reed-group standard triad (Quinn 2018 p. 4641, Quinn 2019
   p. 5965). MOEA-FIND §2.2 uses exactly this triad.
8. Hadka and Reed 2015 Table 2 (p. 361) is the standard format for
   reporting per-objective epsilon precisions. Adopted in §2.7.
9. Reed 2013 p. 441 is the correct citation for epsilon-dominance
   archiving alongside Laumanns et al. 2002, and is cited in §2.3 and
   §2.4 of the manuscript.

Citation corrections applied in this update:

- Kirsch 2013 journal is *Journal of Water Resources Planning and
  Management* 139(4) 396 to 406 (verified, restored).
- Reed 2013 page range is 438 to 456 (verified).
- Herman 2016 journal is *Journal of Water Resources Planning and
  Management* 142(11) article 04016050 (added).

---

## 9. Robustness Frameworks and DMDU Foundations

This section captures foundational robustness-framework and decision-making-under-deep-uncertainty (DMDU) papers held in `manuscript/literature/` that previously lacked dedicated literature-review entries. Together with §1 (search-based generation), §1B (scenario-neutral lineage), §6 (scenario discovery / robustness), and §8 (extraction round 2026-04-14), this completes the lit-review coverage of the literature folder. These entries are intentionally tighter than §1B since the framing-anchor work (`framing_anchor.md`) treats Herman 2015 and Moallemi 2020 in depth already.

### Herman et al. (2015) — Robustness Taxonomy
**Citation:** Herman, J.D., Reed, P.M., Zeff, H.B., and Characklis, G.W. (2015). How Should Robustness Be Defined for Water Systems Planning under Change? *Journal of Water Resources Planning and Management*, 141(10), 04015012. doi:10.1061/(ASCE)WR.1943-5452.0000509

**Foundational robustness-framework taxonomy paper.** Proposes a four-axis taxonomy (Fig. 1) for comparing robustness frameworks:
1. **Alternative generation** — prespecified vs. searched; single- vs. multi-objective
2. **Sampling of states of the world** — key factors assumed vs. discovered; prespecified scenarios vs. design of experiments; outward-from-expected vs. global sampling
3. **Quantification of robustness measures** — expected value, satisficing, regret; domain criterion / uncertainty horizon; univariate vs. multivariate thresholds; deviation from best vs. baseline
4. **Sensitivity analysis controls** — factor mapping (PRIM), factor prioritization (Sobol), local OAT

Demonstrates on the Research Triangle case study (urban water supply, North Carolina) that methodological choices in this taxonomy lead to substantially different planning recommendations. Three load-bearing recommendations:
1. "Decision alternatives should be searched rather than prespecified"
2. "Dominant uncertainties should be discovered through sensitivity analysis rather than assumed"
3. "A carefully elicited multivariate satisficing measure of robustness allows stakeholders to achieve their problem-specific performance requirements"

Introduces or formalizes terms: **a posteriori decision support**, **generate-first-choose-later (GFCL)**, **states of the world (SOWs)**, **domain criterion satisficing measure**, **factor mapping**, **factor prioritization**, **multivariate satisficing**.

**Relevance:** Foundational for vocabulary. Cite for: the four-axis robustness taxonomy (Axis II = "States of the world" is where MOEA-FIND lives), "states of the world / SOWs," "a posteriori decision support," "generate-first-choose-later," "domain criterion satisficing," "factor mapping" / "factor prioritization." See `framing_anchor.md` §1 for the full extraction of Herman 2015 vocabulary applied to the MOEA-FIND positioning.

**Note:** Herman 2015 contemplates only two sub-modes in Axis II ("Prespecified Scenarios" and "Design of experiments") for sampling SOWs. MOEA-FIND argues that *directed search in hazard outcome space* is a third mode not contemplated by this taxonomy (see `framing_anchor.md` §7 paradigm-(iv)).

### Maier et al. (2016) — DMDU Vocabulary Synthesis
**Citation:** Maier, H.R., Guillaume, J.H.A., van Delden, H., Riddell, G.A., Haasnoot, M., and Kwakkel, J.H. (2016). An uncertain future, deep uncertainty, scenarios, robustness and adaptation: How do they fit together? *Environmental Modelling & Software*, 81, 154-164. doi:10.1016/j.envsoft.2016.03.014

Multidisciplinary synthesis paper that articulates how five concepts fit together: **uncertain future**, **deep uncertainty**, **scenarios**, **robustness**, and **adaptation**. Argues these have generally been considered in isolation; this paper provides the connective tissue.

Three complementary paradigms for modeling the future (§2, Fig. 1):
1. **Use of best available knowledge** — single deterministic estimate of the future ("clear enough" future)
2. **Quantification of future uncertainty** — probability distributions over inputs and parameters (aleatory or "Knightian-quantifiable" uncertainty)
3. **Exploration of multiple plausible futures** — set of scenarios with no associated probability or even ranking (deep / Knightian uncertainty)

Important typology distinctions:
- **Aleatory (ontic) uncertainty** — intrinsic uncertainty of natural variability
- **Epistemic uncertainty** — lack of knowledge, or ambiguity (multiple frames of reference)
- **Deep uncertainty (Knightian)** — multiple plausible futures with no probability ranking

Classification of scenario types (§4):
- **Predictive scenarios** — forecast-style, typically with associated probabilities
- **Exploratory scenarios** — descriptive of plausible futures (e.g., GCM scenarios)
- **Normative scenarios** — back-cast or end-state-defined (e.g., visions, target conditions)

**Relevance:** Foundational for vocabulary. Cite for: the three-paradigm classification of future modeling, the aleatory/epistemic/deep uncertainty distinction, and the predictive/exploratory/normative scenario typology. The Fowler et al. (2024) terminology (Table 1) draws partially on this lineage. Maier 2016 is the closest single paper that situates "deep uncertainty" in vocabulary alongside "scenario," "robustness," and "adaptation."

### McPhail et al. (2018) — Robustness Metrics Comparison
**Citation:** McPhail, C., Maier, H.R., Kwakkel, J.H., Giuliani, M., Castelletti, A., and Westra, S. (2018). Robustness Metrics: How Are They Calculated, When Should They Be Used and Why Do They Give Different Results? *Earth's Future*, 6, 169-191. doi:10.1002/2017EF000649

Provides a unifying calculation framework for robustness metrics. Shows that different metrics lead to different rankings of decision alternatives. Categorizes metric suitability based on three decision-maker preferences:
1. Decision context (absolute performance vs. regret)
2. Risk aversion preference
3. Performance-vs-variance vs. higher-moment focus

Tests the framework on three case studies: Adelaide water supply (Australia), Lake Como (Italy), and the Rhine River (Netherlands).

Introduces or formalizes the **Ψ–s notation** for the dual-loop nested structure: Ψ = {ψ_1, ψ_2, ...} are deeply uncertain outer factors; s_{i,j} ~ G(ψ_i) are inner stochastic realizations conditional on ψ_i. Robustness R(a) = Φ({y}) of decision a is aggregated across the inner ensemble for each outer ψ. This notation is the formal anchor for `framing_anchor.md` §7 paradigm-(iii)/(iv) decomposition.

**Relevance:** Foundational for vocabulary. Cite for: the Ψ–s nested-loop notation, the absolute-vs-regret robustness distinction, and the three-axis decision-maker-preference taxonomy of robustness metrics. The McPhail 2018 framework explains why different aggregation choices lead to different rankings — directly supporting the Hadjimichael 2024 "consequential scenarios depend on aggregation choices" argument cited in MOEA-FIND §4.5.

### Marchau et al. (2019) — Decision Making under Deep Uncertainty (DMDU)
**Citation:** Marchau, V.A.W.J., Walker, W.E., Bloemen, P.J.T.M., and Popper, S.W. (Eds.) (2019). *Decision Making under Deep Uncertainty: From Theory to Practice*. Springer, Cham. doi:10.1007/978-3-030-05252-2

Edited volume that consolidates the DMDU community's vocabulary and methodological framework. Open-access reference text for the field. Chapter coverage:
- Robust Decision Making (RDM) — Lempert chapter
- Dynamic Adaptive Policy Pathways (DAPP) — Haasnoot, Kwakkel, Walker chapter
- Info-Gap Theory — Ben-Haim chapter
- Engineering Options Analysis — de Neufville chapter
- Many-Objective Robust Decision Making (MORDM) — Kasprzyk, Reed, Kirsch chapter
- Decision Scaling — Brown chapter
- Adaptive Delta Management — Bloemen, Hijdra, Kwakkel chapter

**Relevance:** Reference text for vocabulary cross-references. Cite for: definitive expositions of RDM, DAPP, MORDM, and Decision Scaling as named methodological families. Useful for grounding "DMDU" as the umbrella term.

### Herman et al. (2014) — MORDM Multistakeholder Origin
**Citation:** Herman, J.D., Zeff, H.B., Reed, P.M., and Characklis, G.W. (2014). Beyond optimality: Multistakeholder robustness tradeoffs for regional water portfolio planning under deep uncertainty. *Water Resources Research*, 50(10), 7692-7713. doi:10.1002/2014WR015338

Demonstrates that "optimal" Pareto-approximate solutions found under expected future conditions may suffer significant performance degradation under modest deviations in deeply uncertain factors. Introduces the **multistakeholder MORDM framework** that blends many-objective search with uncertainty analysis to discover key tradeoffs between water supply alternatives and their robustness. Demonstrated on the Research Triangle four-utility case study (Raleigh, Durham, Cary, Chapel Hill in North Carolina). Uses PRIM to identify which uncertain factors drive individual and collective vulnerabilities.

**Relevance:** Cite for the multistakeholder MORDM framework as MOEA-FIND's downstream consumer. Demonstrates the explicit pairing of many-objective search with PRIM-based scenario discovery — the Reed-group canonical workflow that MOEA-FIND ensembles feed into.

### Trindade et al. (2017) — MORDM Drought Risk Triggers
**Citation:** Trindade, B.C., Reed, P.M., Herman, J.D., Zeff, H.B., and Characklis, G.W. (2017). Reducing regional drought vulnerabilities and multi-city robustness conflicts using many-objective optimization under deep uncertainty. *Advances in Water Resources*, 104, 195-209. doi:10.1016/j.advwatres.2017.03.023

Advances the multistakeholder MORDM framework by showing that **adaptive risk-of-failure action triggers** must be stressed with a comprehensive sample of deeply uncertain factors *during* the computational search phase (not only after). Search-under-deep-uncertainty fundamentally changes perceived performance tradeoffs and substantially improves robustness compared to search-under-historical-conditions followed by ex-post robustness evaluation. Demonstrated on the Research Triangle case study. Uses MORDM to identify how cooperative water transfers, financial risk-mitigation tools, and coordinated regional demand management must be employed jointly.

**Relevance:** Cite for the methodological advance that uncertainty must enter the search loop, not just the post-hoc evaluation. This argument applies analogously to MOEA-FIND: scenario generation must structure the hazard space *during* generation, not as a post-hoc filtering of a generic ensemble. The Trindade 2017 case study and the related Gold et al. (2023) and Lau et al. (2023) work form a continuous Reed-group methodological lineage that MOEA-FIND ensembles serve.

### Hadjimichael et al. (2023) — Multi-Actor Multi-Impact Scenario Discovery (Preprint)
**Citation:** Hadjimichael, A., Reed, P.M., Quinn, J.D., Vernon, C.R., and Thurber, T.B. (2023). Multi-actor, multi-impact scenario discovery of compound human-natural system futures. ESSOAr preprint (companion to Hadjimichael 2024).

Preprint companion to Hadjimichael 2024 (which is fully covered in §8.1). Same FRNSIC framework. Mentioned here for completeness — the published version is the canonical citation; the preprint is held in `literature/` for the additional methodological detail and figures that did not survive editing for the published version. Defer to Hadjimichael 2024 entry in §8.1 for content.

---
