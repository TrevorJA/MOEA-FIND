# Comprehensive Literature Review: MOEA-FIND

*Date: 2026-04-12*
*Purpose: Ground all experimental design and method choices in published literature*

---

## 1. Search-Based Synthetic Streamflow Generation

### Borgomeo et al. (2015) — "Numerical Rivers"
**Citation:** Borgomeo, E., et al. (2015). Numerical rivers: A synthetic streamflow generator for water resources vulnerability assessments. Water Resources Research, 51(7).

Pioneered treating streamflow generation as a combinatorial optimization problem. Uses simulated annealing to shuffle values in initialized synthetic time series until sequences match user-defined objectives (target statistical properties). Enables representation of climate-induced changes in persistence, variability, and extremes. Single-objective formulation: one target, one trace per optimization run.

**Relevance:** Direct intellectual predecessor. MOEA-FIND generalizes to multi-objective, where the Pareto front is the ensemble.

### Zaniolo et al. (2023) — FIND
**Citation:** Zaniolo, M., Fletcher, S., and Mauter, M. (2023). FIND: A synthetic weather generator to control drought Frequency, INtensity, and Duration. Environmental Modelling and Software.

Extends Borgomeo's approach with direct, independent control of three drought properties (FID) via Standardized Streamflow Index. Uses simulated annealing. Designed for bottom-up vulnerability analysis: "train and test drought strategies under historical and plausible future drought conditions." MATLAB implementation available.

**Relevance:** Closest precedent. MOEA-FIND addresses FIND's single-objective limitation and automates ensemble generation.

### Wheeler et al. (2024) — Multi-Site Extension
**Citation:** Wheeler, K.G., et al. (2024). Multisite Nonparametric Stochastic Streamflow Generation for the Eastern Nile Basin. Journal of Hydrologic Engineering, 30(1).

Extends simulated annealing generation to multiple sites using cross-site correlation matrix in the objective function. Includes Hurst coefficient for long-range dependence. Demonstrates scalability to complex transboundary basins.

**Relevance:** Shows search-based generation works at multi-site scale. MOEA-FIND's block bootstrap with shared indices is a simpler approach to spatial correlation.

### Herman et al. and Reed Group — MORDM Applications
The Reed group (Cornell) extensively employs synthetic streamflow generation within the Many-Objective Robust Decision Making (MORDM) framework. Synthetic generation enables: (1) preserving effects of deeply uncertain factors, (2) evaluating policies across plausible futures, (3) discovering system vulnerabilities via ensemble simulation. OpenMORDM software documents this workflow.

**Relevance:** MOEA-FIND fits naturally as the scenario generation component of MORDM, providing structured drought coverage that unstructured random ensembles lack.

---

## 2. Kirsch-Nowak Block Bootstrap Methods

### Kirsch et al. (2013)
**Citation:** Kirsch, B.R., et al. (2013). Evaluating the impact of alternative hydro-climate scenarios on transfer agreements: Practical improvement for generating synthetic streamflows. Journal of Water Resources Planning and Management, 139(4), 396-406.

The standard block bootstrap for synthetic streamflow. Operates in log-space: (1) log-standardize historical flows, (2) moving block bootstrap in log-space, (3) apply Cholesky decomposition for historical correlation, (4) de-standardize and back-transform. Addresses Fractional Gaussian Noise limitations. Widely implemented.

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
- Zaniolo et al. (2023) — FIND comparison
- Svensson et al. (2017) — Drought threshold/metrics

### Parametric CDF (Experiments 3.1-3.3)
- Svensson et al. (2017) — Kappa distribution
- PAR(p)-vine copula (Brechmann et al. 2017)
- Papalexiou and Koutsoyiannis — Phase randomization
- Mixed copula-KDE (2025)

### Multi-Site DRB (Experiments 4.1-4.3)
- Wheeler et al. (2024) — Multi-site generation
- Nowak et al. (2010) — Multi-site disaggregation
- Hadjimichael et al. (2020) — Stakeholder robustness
- BART (Chipman et al. 2010; Water Programming 2024)
