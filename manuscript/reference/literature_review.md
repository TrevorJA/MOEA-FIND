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

### Herman et al. and Reed Group — MORDM Applications
The Reed group (Cornell) extensively employs synthetic streamflow generation within the Many-Objective Robust Decision Making (MORDM) framework. Synthetic generation enables: (1) preserving effects of deeply uncertain factors, (2) evaluating policies across plausible futures, (3) discovering system vulnerabilities via ensemble simulation. OpenMORDM software documents this workflow.

**Relevance:** MOEA-FIND fits naturally as the scenario generation component of MORDM, providing structured drought coverage that unstructured random ensembles lack.

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
