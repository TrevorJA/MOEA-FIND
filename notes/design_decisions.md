# Design Decisions and Open Questions

*Living document. Updated as decisions are made.*

---

## DD-01: Trace Length and Drought Event Design

### The question

Should MOEA-FIND generate long traces (e.g., 78-year records matching historical length) or short traces containing only 1-2 drought events?

### Options

**Option A: Long traces (78 years)**
- Matches the simulation period of Pywr-DRB
- Drought metrics (frequency, mean duration, mean intensity) are computed as averages over the full trace
- Each trace contains multiple drought events whose properties are summarized
- Decision variables: N_years * 12 = 936 DVs for 78 years
- Pro: Directly usable in Pywr-DRB without modification
- Con: High dimensionality. "Mean drought severity" averages over multiple events, obscuring individual event characteristics. Borg must search a 936-dim space.

**Option B: Short traces (5-15 years, containing 1-2 droughts)**
- Each trace is a focused drought scenario
- Drought metrics describe the specific event(s), not long-term averages
- Decision variables: N_years * 12 = 60-180 DVs (much more tractable)
- Pro: Much fewer DVs. Drought characteristics directly describe the event, not an average. Faster Borg convergence. More interpretable: "this trace IS a specific drought."
- Con: Cannot be directly fed into a 78-year Pywr-DRB simulation. Would need to embed the short drought trace within a longer context (e.g., sandwich between normal-flow years). Initial conditions matter more for short traces.

**Option C: Hybrid (moderate length, 20-30 years)**
- Contains 2-4 drought events depending on frequency target
- Enough length for meaningful frequency statistics
- DVs: 240-360 (reasonable for Borg)
- Compromise between A and B

### Considerations

- For the DRB case study, Pywr-DRB expects a specific simulation length and needs spin-up time for reservoir storage to equilibrate. Short traces require careful handling of initial conditions.
- The "generate short events, embed in context" approach has precedent in the climate stress testing literature (e.g., splicing drought events into baseline sequences).
- If traces are short and contain only 1 drought, "frequency" as an objective becomes meaningless. The objectives would shift to event-level metrics: duration, peak intensity, cumulative severity, onset timing.
- Short traces are more natural for "stress testing" (what's the worst plausible drought?). Long traces are more natural for "robustness evaluation" (how does the policy perform across varied conditions?). These are different questions.

### Decision (2026-04-12)

**Both framings should be explored.** The publication should demonstrate both levels of temporal application:

1. **Event-level traces (1-8 year events):** More novel contribution. Event-specific metrics (duration, peak intensity, cumulative severity) avoid hiding insights behind aggregate averages. DVs are tractable. This is the more interesting and novel framing. Closest precedent is the climate stress testing literature on splicing designed drought events into baseline sequences. SDF (Severity-Duration-Frequency) curves provide the established framework for joint event-level characterization.

2. **Moderate-length traces (~30 years):** Aggregate metrics (mean duration, frequency, mean intensity) are better supported by existing literature. Constraining to ~30 years keeps DVs manageable (~360 monthly or ~120 seasonal-block DVs). Needed for literature comparability. Matches the framing in Borgomeo et al. (2015), Zaniolo et al. (2023, FIND), and Wheeler et al. (2024).

For the proof-of-concept paper, demonstrate both and compare what each reveals. The event-level framing is the stronger contribution; the aggregate framing provides validation against established approaches.

---

## DD-02: Decision Variable Formulation

### Options

**Option A: Discrete bootstrap indices**
- Each DV is an integer index selecting a historical month
- Search space: integers in [0, N_months_historical)
- Pro: Direct, no mapping artifacts
- Con: Discrete space. Borg's continuous operators (SBX, DE) produce real values that must be rounded. May reduce search efficiency.

**Option B: CDF probabilities**
- Each DV is p in [0, 1], mapped to a historical month via inverse empirical CDF
- Search space: continuous [0, 1]^d
- Pro: Continuous, natural for Borg. Small perturbations in p produce small perturbations in selected month (mostly).
- Con: Non-uniform mapping. Months near the median have higher density in CDF space. Quantization at CDF boundaries.

**Option C: Normalized index (continuous relaxation)**
- Each DV is x in [0, 1], mapped to month index via floor(x * N_months_historical)
- Search space: continuous [0, 1]^d
- Pro: Simple, uniform mapping
- Con: Small perturbation in x can cause a jump if crossing a boundary

**Option D: Parametric CDF sampling (extends beyond historical)**
- Fit parametric distributions (gamma, log-normal, kappa) to historical monthly flows
- Each DV is p in [0, 1], mapped to a flow value via the inverse fitted CDF
- Search space: continuous [0, 1]^d
- Pro: Generates flows outside the historical range. Smooths the empirical distribution. Natural for Borg.
- Con: Requires defensible distribution fitting. Temporal/spatial correlation must be imposed separately (via copula or Cholesky). Adds complexity.
- **Literature support:** Kappa distribution (Svensson et al. 2017) is flexible with bounded support below zero. Vine copulas (D-vine, periodic) handle temporal lag-1/lag-2 and spatial dependence (PAR(p)-vine models, 2017-2023). This is current state-of-the-art for synthetic streamflow generation.

**Option E: KDE-smoothed CDF**
- Estimate kernel density for each month's historical flows
- Each DV is p in [0, 1], mapped via inverse KDE CDF
- Pro: Nonparametric but smoother than bootstrap. Extends support via kernel bandwidth.
- Con: Less extrapolation than parametric. Bandwidth selection affects tail behavior.
- **Literature support:** Mixed copula-KDE framework (2025, Stoch. Environ. Res. Risk Assess.)

**Option F: Bootstrapped CDF uncertainty**
- Bootstrap the historical data to generate multiple empirical CDF estimates with confidence bands
- DVs map to quantiles of the ensemble of CDFs, not a single CDF
- Pro: Formally accounts for estimation uncertainty. Ensemble of CDFs naturally extends support.
- Con: Adds a layer of complexity. Need to decide how to sample across the CDF ensemble.
- **Literature support:** BEUM approach (2022, Stoch. Environ. Res. Risk Assess.)

### Considerations

- Borgomeo (2015) and Wheeler (2024) used simulated annealing with discrete swap operations (swap one month for another), which avoids the continuous/discrete issue
- Borg is designed for continuous variables. Option B, D, or E preferred.
- With block bootstrapping, DVs may represent block start indices rather than individual months, reducing dimensionality further
- **Critical issue (from user feedback):** Pure bootstrap resampling limits generated droughts to recombinations of historical patterns. Options D-F allow generating droughts more extreme than anything in the record. This is a key methodological question for the paper.
- **Phased approach:** Start with bootstrap (Options B/C) for POC. Then implement parametric CDF (Option D) with vine copula for correlation. Compare the two in terms of drought space coverage, especially in the tails.
- **Kappa distribution** (Svensson et al. 2017) is the most defensible parametric choice for monthly streamflow: bounded below at zero, flexible tail behavior, widely used.
- **Vine copulas** (D-vine with lag-1/lag-2) are state-of-the-art for temporal dependence when using parametric marginals (PAR(p)-vine, 2017-2023 literature)
- **Phase randomization** (Papalexiou and Koutsoyiannis 2019) is an alternative that preserves spectral properties while allowing parametric tail extension

### Decision (2026-04-12)

**Two-track approach.** Bootstrap resampling (Options B/C) for initial implementation and POC. Parametric CDF with vine copula (Option D) as the advanced formulation for the publication, enabling generation of droughts beyond the historical envelope. Both tracks use Borg with the same Manhattan norm trick. The comparison between them is itself a contribution: does parametric sampling meaningfully expand the achievable drought space?

**Literature grounding:**
- Bootstrap track: Kirsch et al. (2013) block bootstrap; Borgomeo et al. (2015) search-based generation with swap operations; Wheeler et al. (2024) multi-site with discrete swaps.
- Parametric track: Svensson et al. (2017) kappa distribution for monthly marginals; Brechmann et al. (2017) PAR(p)-vine copula for temporal dependence; Li et al. (2023) mixed D-vine conditional quantile model; periodic vine copulas (2019) for seasonal variation.
- Middle-ground options: Mixed copula-KDE framework (2025, Stoch. Environ. Res. Risk Assess.) for nonparametric smoothing; BEUM bootstrapped CDF uncertainty (2022); Papalexiou (2019) phase randomization preserving spectral properties.

---

## DD-03: Bootstrap Granularity (Monthly vs. Block-Level DVs)

### The question

Should each DV select an individual month, or should DVs select blocks of months?

### Options

**Option A: Monthly DVs**
- One DV per month of the synthetic trace
- Maximum flexibility in recombination
- Destroys temporal autocorrelation within months (unless constrained)

**Option B: Block-level DVs**
- One DV selects a starting month; the next B months are copied as a block
- DVs: N_years * 12 / B (e.g., with B=3, 936/3 = 312 DVs for 78-year trace)
- Preserves within-block temporal structure automatically
- Reduces dimensionality
- Block boundary artifacts (discontinuities)

**Option C: Seasonal blocks**
- Fixed block structure aligned to seasons (e.g., DJF, MAM, JJA, SON)
- DVs: N_years * 4 (e.g., 312 for 78 years, or 40-60 for 10-15 year traces)
- Preserves seasonal autocorrelation perfectly within blocks
- Natural for hydrological applications where seasons are the fundamental unit

### Considerations

- The Kirsch generator uses block bootstrapping with variable block size. Standard is B=1-3 months.
- For drought generation, preserving multi-month persistence is critical. Longer blocks help.
- For short traces (DD-01 Option B), seasonal blocks (Option C) would give very few DVs (e.g., 10-year trace = 40 seasonal DVs). Highly tractable for Borg.
- Block-level DVs with blocks of years (B=12) would reduce a 78-year trace to 78 DVs. Each DV selects which historical year to resample. This preserves annual patterns perfectly but limits recombination to annual shuffling.

### Decision (RESOLVED 2026-04-13)

**Kirsch pipeline uses B=1 (monthly) + Cholesky, not mechanical blocks.** The Kirsch method (Kirsch et al. 2013) uses independent monthly resampling (B=1) and imposes temporal correlation post-hoc via Cholesky decomposition of the historical monthly correlation matrix in normal-score space. The Dec-Jan boundary is handled separately via a shifted Y_prime matrix. This is distinct from a general block bootstrap where B>1 preserves within-block autocorrelation mechanically.

For MOEA-FIND, the KirschBorgWrapper (see DD-06) delegates to SynHydro's validated Kirsch implementation, which includes the full pipeline: normal-score transform, Cholesky factor application, inverse transform, and Dec-Jan handling. This approach avoids reimplementing these critical steps.

Alternative approaches exist: larger blocks (B=3-6) can preserve within-block temporal structure mechanically without Cholesky, producing simpler but less flexible implementations. However, mechanical blocks are not the Kirsch method and introduce block-boundary discontinuities. The full Kirsch pipeline (B=1 + Cholesky) is preferred for publication-quality work.

**Literature on block size selection:** HHJ (Hall, Horowitz, Jing) cross-validation and NPPI (Nonparametric Plug-In) methods estimate optimal block length for bootstrap variance reduction. These apply to mechanical block bootstrap, not Kirsch's Cholesky pipeline.

---

## DD-04: Objective Set Design

### The question

Which drought characteristics should be objectives, and how should they be defined?

### Options for drought characteristics (each would be an objective + Manhattan norm)

**Event-level metrics (for short traces with 1-2 droughts):**
- Duration: number of consecutive months below threshold
- Peak intensity: maximum single-month deficit below threshold
- Cumulative severity: total deficit below threshold over the event
- Onset timing: month of year when drought begins (seasonality)
- Recovery rate: how quickly flow returns to normal after the event

**Trace-level metrics (for long traces with multiple droughts):**
- Mean/max drought duration across all events
- Mean/max average severity (total deficit divided by duration)
- Drought frequency (events per decade)
- Inter-event time (mean time between droughts)

### Threshold definition

- Fixed percentile of historical: P20, P30?
- Basin-specific threshold (e.g., FFMP drought trigger levels)?
- Variable threshold by season?

### Anti-ideal point placement

For the Manhattan norm objective, D* should be at the "worst plausible" corner:
- Maximum frequency, longest duration, highest intensity
- These maxima should be based on some multiple of historical maxima (e.g., 2x) or physical limits
- If too far from the actual Pareto front, epsilon dominance becomes less effective
- If too close, the Pareto front may not span enough of the drought space

### Decision

**TBD.** Depends critically on DD-01 (trace length determines whether we use event-level or trace-level metrics).

**Literature on drought metrics and thresholds:**
- SSI (Standardized Streamflow Index) used by Zaniolo et al. (2023, FIND) for drought control, with onset at SSI <= -0.84 (20th percentile).
- P20 monthly flow (20th percentile of flow duration curve) is a defensible default fixed threshold. Variable daily threshold (VDT, 95th percentile via 31-day moving window) accounts for seasonality but adds complexity.
- SDF (Severity-Duration-Frequency) curves from bivariate copula-based frequency analysis inform the range of characteristics MOEA-FIND should target.
- Multi-dimensional drought characterization integrates meteorological (weeks), agricultural (months), and hydrological (months-years) timescales. MOEA-FIND focuses on hydrological drought.
- SPEI (Vicente-Serrano et al.) may be preferable to SPI for warming scenarios as it includes temperature via reference ET.

---

## DD-05: Constraint Design

### The question

What constraints should enforce trace plausibility, and how tight should they be?

### Candidates

1. **Lag-1 autocorrelation:** |rho_1(Q_S) - rho_1(Q_H)| < delta
2. **Annual statistics:** |mean(Q_S) - mean(Q_H)| < delta (non-drought years)
3. **Seasonal cycle:** monthly means within X% of historical
4. **Hurst coefficient:** long-range dependence (from Wheeler 2024)
5. **Cross-site correlation:** for multi-site generation
6. **Non-negative flows:** Q_S >= 0 everywhere
7. **No impossible extremes:** max(Q_S) < some physical limit

### Considerations

- Too many constraints: Borg's constraint handling reduces the feasible population. May prevent finding good coverage of drought space.
- Too few constraints: Degenerate traces (all zeros, unrealistic patterns)
- The block bootstrap inherently preserves much temporal structure. If blocks are large enough (seasonal), constraints 1-3 may be automatically satisfied.
- Constraint violations can be evaluated cheaply alongside drought metrics
- Starting with minimal constraints and adding if needed is more robust than starting strict

### Decision

**TBD.** Start with minimal constraints (non-negative flows, basic autocorrelation). Add incrementally based on plausibility assessment of generated traces.

**Literature on constraints:** Wheeler et al. (2024) includes Hurst coefficient for long-range dependence and cross-site correlation matrix in the objective function (not as a constraint, but serves the same purpose). Borgomeo et al. (2015) used autocorrelation matching as part of the simulated annealing objective. The block bootstrap (Kirsch et al. 2013) inherently preserves within-block temporal structure and, with shared indices, cross-site correlation. Larger blocks reduce the need for explicit autocorrelation constraints. The Cholesky decomposition step in the Kirsch pipeline further enforces historical correlation structure.

---

## DD-06: SynHydro Integration vs. Standalone Implementation

### The question

Should MOEA-FIND call SynHydro's internals, or implement its own bootstrap generator?

### Options

**Option A: SynHydro integration**
- Import SynHydro's bootstrap/resampling functions
- Need SynHydro to expose a lower-level API (accept external indices, return trace)
- Pro: reuse validated code, consistent with the SynHydro publication
- Con: API coupling. SynHydro may not expose the right interface. Changes to SynHydro affect MOEA-FIND.

**Option B: Standalone implementation**
- Implement block bootstrap in MOEA-FIND's own `generator.py`
- Simple: it's just resampling with shared indices across sites
- Pro: self-contained, no external dependency, can optimize for MOEA-FIND's specific needs (e.g., CDF-probability interface)
- Con: code duplication. Must validate independently.

**Option C: Thin wrapper**
- MOEA-FIND defines the bootstrap logic but uses SynHydro for data loading, preprocessing, and daily disaggregation
- Separation of concerns

### Considerations

- The block bootstrap itself is simple (<100 lines). The value of SynHydro is in the full pipeline (data loading, deseasonalization, Cholesky correlation, daily disaggregation).
- For a publication-quality implementation, Option C (wrapper using SynHydro for preprocessing) seems best.
- For rapid prototyping, Option B gets us moving fastest.

### Decision (RESOLVED 2026-04-13)

**Option C (KirschBorgWrapper using SynHydro's validated pipeline).** Analysis revealed that a standalone `BlockBootstrapGenerator` was missing three critical components of the Kirsch pipeline: normal-score transform, Cholesky decomposition, and Dec-Jan boundary handling. These components are essential for preserving historical autocorrelation and plausibility. Rather than reimplementing them independently, the wrapper integrates with SynHydro's validated Kirsch generator.

The wrapper supports **two DV injection modes**, corresponding to the two sampling formulations (see DD-02):

**Mode 1: Index injection (discrete bootstrap)**
Borg DVs are continuous [0,1] values mapped to discrete historical year indices via `floor(p * n_years_hist)`. These indices replace SynHydro's random year selection. The full Kirsch pipeline (deseasonalization, normal-score transform, Cholesky, Dec-Jan handling) runs on the selected historical data. This preserves the block-resampling character of Kirsch while allowing Borg to steer which years are sampled.

**Mode 2: Residual injection (CDF-based)**
Borg DVs are CDF probabilities that map to standardized residuals via the empirical or KDE-smoothed CDF. These residuals are then fed into SynHydro's Cholesky + normal-score pipeline to impose historical correlation structure. This keeps the continuous DV landscape smoother for Borg's operators while gaining the correlation structure that a naive standalone generator would lack.

Both modes use the same SynHydro backend for correlation imposition. The choice between them is an empirical question for the paper (does one produce better Borg convergence or coverage quality?).

**SynHydro API requirement:** Expose the post-resampling pipeline (normal-score transform, Cholesky multiplication, inverse normal-score, Dec-Jan handling) as a separate public method that accepts externally provided residuals or year indices. Small, clean addition to the existing `KirschGenerator` class.

**Wrapper interface sketch:**

```python
class KirschBorgWrapper:
    """Wraps SynHydro KirschGenerator for Borg DV injection.
    
    Args:
        kirsch: Fitted SynHydro KirschGenerator instance.
        mode: "index" or "residual".
        n_months_out: Number of output months (trace length * 12).
    """
    
    def generate(self, dvs: np.ndarray) -> np.ndarray:
        """Generate a synthetic trace from Borg decision variables.
        
        Args:
            dvs: Array of [0, 1] values, length n_months_out.
                In "index" mode: mapped to year indices.
                In "residual" mode: mapped to standardized residuals.
        
        Returns:
            Monthly synthetic flows, shape (n_years_out, 12).
        """
        if self.mode == "index":
            indices = self._dvs_to_year_indices(dvs)
            return self.kirsch.generate_from_indices(indices)
        else:
            residuals = self._dvs_to_residuals(dvs)
            return self.kirsch.generate_from_residuals(residuals)
```

**Literature context:** Borgomeo et al. (2015) and Wheeler et al. (2024) both implemented standalone generators with simulated annealing. The Kirsch bootstrap itself is simple (<100 lines), so standalone implementation is low-risk. The value of SynHydro integration is in the full pipeline: data loading, log-space deseasonalization, Cholesky correlation decomposition, and daily disaggregation (Nowak et al. 2010). For the MOEA-FIND publication, the monthly generation mechanism needs to be clearly documented. The KirschBorgWrapper provides this clarity while eliminating code duplication and ensuring correctness of the Cholesky + normal-score pipeline.

---

## DD-07: Borg Interface

### The question

How to interface with Borg MOEA?

### Options

**Option A: borg.py (Python wrapper)**
- Standard Borg Python API used in NYCOptimization
- Requires compiled Borg C library
- Single-master, master-worker, or multi-master modes available

**Option B: platypus (Python MOEA library)**
- Open-source, includes Borg-like algorithms (epsilon-MOEA, NSGA-II)
- No licensing issues
- May not have identical epsilon-dominance behavior to Borg

**Option C: pymoo**
- Modern Python MOEA library
- NSGA-III, MOEA/D, etc.
- No epsilon-dominance archiving equivalent to Borg

### Considerations

- The Manhattan norm trick relies on epsilon-dominance for uniformity. Borg's epsilon-box archiving is essential to the method.
- platypus has an epsilon-MOEA implementation that may work similarly. Worth testing.
- For the publication, using Borg (the standard in the water resources community) is strongest.
- For the proof-of-concept, platypus may be sufficient and avoids Borg licensing.

### Decision

**TBD.** Use platypus for POC experiments, Borg for publication results. Verify that platypus epsilon-MOEA produces similar uniformity.

**Literature context:** Hadka and Reed (2013) demonstrated Borg's superiority on 33 benchmark instances (DTLZ, WFG, CEC 2009). Three key mechanisms make Borg essential for MOEA-FIND: (1) epsilon-dominance archiving produces near-uniform tiling of the Manhattan-norm hyperplane, (2) epsilon-progress restart triggering escapes premature convergence in the high-dimensional DV space, (3) adaptive operator selection among six variational operators (SBX, DE, PCX, UNDX, SPX, PM, UM) adapts to the problem landscape. Per user direction, Borg is the required optimizer for this method; generalization to other MOEAs is out of scope. NSGA-III (Deb and Jain 2014) and MOEA/D (Zhang and Li 2007) achieve coverage through different mechanisms (reference points, decomposition) that lack the specific hyperplane-tiling property.

---

## DD-08: Zero-Drought Solutions on the Pareto Front

### Observation (2026-04-13)

The Kirsch (residual) mode produces exactly 1 Pareto solution at (0, 0) in drought space (mean_duration=0, mean_avg_severity=0), representing a trace with no detected SSI-3 drought events.

### Analysis

This is mathematically correct and expected. The solution sits at the corner of the Manhattan-norm hyperplane: J1=0, J2=0, J3=sum(anti_ideal)=64.55. It is non-dominated because no other solution can achieve J1=0 AND J2=0 simultaneously. The Pareto front correctly spans from the "no drought" corner to the extreme drought corner, providing full coverage of the feasible drought characteristic space.

The (0, 0) point is a natural endpoint of the ensemble. Whether it is useful depends on the downstream application. For robustness evaluation, a "no drought baseline" scenario may actually be informative.

### Options if exclusion is desired

**Option A: Constraint (n_events >= 1).** Add a plausibility constraint requiring at least one drought event in the generated trace. Removes the zero-drought corner while keeping the rest of the Pareto front intact. Simple and principled.

**Option B: Shift objective bounds.** Set minimum drought metric values above zero (e.g., duration >= 1 month). This shifts the ideal point away from (0,0). More complex, affects epsilon grid.

**Option C: Do nothing.** At higher NFE budgets with more Pareto solutions, a single corner point is negligible. The ensemble is not degraded by including one benign scenario.

### Decision

**Deferred.** Not a priority at current stage. Revisit if zero-event solutions become numerous at higher NFE or if downstream analysis is affected. Option A (n_events >= 1 constraint) is the recommended fix if needed.

---

## DD-09: Coverage Baseline Comparison Bounds

### Observation (2026-04-13)

The L2-star discrepancy comparison between MOEA-FIND and LHS/Sobol baselines (Fig 3 in diagnostics) is not apples-to-apples. LHS and Sobol baselines fill the full bounding box [0, anti_ideal] in drought space, while MOEA-FIND only covers the feasible region (a subset determined by the generator's capabilities and the physics of drought).

### Issue

MOEA-FIND's L2* (0.47 index, 0.28 residual) appears much worse than LHS/Sobol (~0.03), but this comparison is misleading. LHS/Sobol sample the full rectangle, including drought configurations that may be physically infeasible (e.g., very long duration with very low severity). MOEA-FIND only places solutions where the generator can actually produce droughts.

### Options

**Option A: Constrained baselines.** Generate LHS/Sobol samples in the full rectangle, then filter to only those that fall within the convex hull (or some envelope) of the MOEA-FIND Pareto front. Compute discrepancy on the shared feasible region.

**Option B: Empirical CDF transform.** Transform both MOEA-FIND and baseline points to the empirical marginal CDF of the MOEA-FIND Pareto front, then compute discrepancy. This normalizes for the feasible region shape.

**Option C: Relative discrepancy.** Report MOEA-FIND discrepancy relative to the LHS/Sobol discrepancy at the same sample size, acknowledging the different domains. Focus the narrative on nearest-neighbor CV (which is domain-agnostic for spacing regularity) rather than L2*.

**Option D: Domain-specific baselines.** Generate a large Kirsch library (10K+ random traces), compute their drought characteristics, then subsample using LHS in drought characteristic space. This is the most directly comparable baseline (same feasible region, same generator).

### Decision (updated 2026-04-13)

**Option D implemented** as `src/library.py` with experiment script `experiments/kirsch_ensemble/run_library_baseline.py`. Supports both LHS and Sobol subsampling in drought characteristic space. Subsampling uses nearest-neighbor matching from design points to library members.

Infrastructure is ready. Large-scale library generation (10K+ traces) deferred to HPC or dedicated local run.

### Literature Context

**Synthetic streamflow generation in water resources:**
Kirsch et al. (2013) introduced the block bootstrap stochastic streamflow generator with Cholesky correlation preservation (B=1 monthly resampling + normal-score transform). The Kirsch-Nowak generator is widely used in Reed group MORDM studies to produce synthetic streamflow ensembles for robustness evaluation. Typical practice generates ensembles of 100-10,000 traces and uses them directly for policy re-evaluation without structured subsampling in drought characteristic space.

**How existing studies actually use synthetic ensembles:**
Hadjimichael et al. (2020) samples States of the World via LHS in *parameter space* (HMM transition probabilities, mean/variance multipliers), then generates synthetic realizations per SOW. The LHS operates on generator parameters, not on drought characteristics post-hoc. Quinn et al. (2018) uses CMIP5 climate projections with bias-corrected statistical downscaling, not synthetic library subsampling. Herman et al. (2016) modifies stochastic generation to control drought properties during generation via SSI, rather than selecting post-hoc.

**Direct optimization approaches:**
Borgomeo et al. (2015) uses simulated annealing to shuffle synthetic sequences until they match target statistical properties (single-objective). Zaniolo et al. (2024) extends this with block bootstrap + direct SSI-based drought control, allowing independent specification of frequency, intensity, and duration. Both generate on-demand to match targets rather than subsampling from pre-existing libraries.

**Space-filling experimental design and scenario subsampling:**
LHS (McKay et al., 1979) and Sobol quasi-random sequences (Sobol, 1967) are standard space-filling designs widely used for sampling uncertainty factors in sensitivity analysis and MORDM (reviewed in Pianosi et al., 2016). Bonham et al. (2024, "Subsampling and space-filling metrics to test ensemble size for robustness analysis," Environ. Model. Softw., 172, 105933) is the closest published precedent: they use conditioned LHS (cLHS) to subsample from a pre-generated 500-scenario library in uncertainty input space (streamflow features, demand, initial storage) for the Colorado River Basin. They evaluate subset quality via space-filling metrics (MST mean edge length, minimum distance) and show that MSTmean predicts robustness ranking accuracy (R2 = 0.77-0.91). However, their subsampling operates in the uncertainty *input* space (generator parameters and streamflow statistics), not in drought *outcome* space (drought duration, severity, frequency). Their focus is on computational efficiency (how few scenarios preserve policy rankings), not on structured coverage of drought characteristics for vulnerability discovery.

**Gap that MOEA-FIND addresses:**
The specific workflow of subsampling a synthetic streamflow library in drought *characteristic* space (duration, severity, frequency) has not been published. Bonham et al. (2024) subsample in input space; Hadjimichael et al. (2020) sample in parameter space. Neither targets structured coverage of drought outcomes. MOEA-FIND goes further by directly optimizing trace generation for coverage via the Manhattan norm trick, bypassing the library entirely. The comparison between library-subsample in drought space (post-hoc selection, implemented here as a novel baseline) and MOEA-FIND (direct optimization for coverage) tests whether optimization during generation outperforms post-hoc selection from a finite random library. This comparison is a core contribution of the paper (RQ4).

---

## DD-10: Epsilon-Dominance Coverage Quality in Higher Dimensions

### Observation (2026-04-13, Session 14)

The 3D analytic validation test (Experiment 1.2 complete) reveals that epsilon-dominance tiling on the (k-1)-simplex produces moderate, not excellent, coverage uniformity in the full DV space when k=3.

**Results:**
- 1362 Pareto-optimal solutions on the 2-simplex (triangular surface in 3D objective space)
- Hyperplane constraint verified: solutions satisfy sum(J_i) = constant to machine precision (~10^-16 std)
- Full range coverage: all three objectives span [-3, 3]
- Nearest-neighbor coefficient of variation (NN_CV): 0.42 for Pareto front
- L2-star discrepancy (L2*): 0.038 for Pareto front

**Comparison to baselines on the same [0, anti_ideal]^3 bounding box:**
- LHS (1362 samples): NN_CV = 0.37, L2* = 0.006
- Sobol (1362 samples): NN_CV = 0.28, L2* = 0.001
- Random (1362 samples): NN_CV = 0.50, L2* = 0.042

**Interpretation:** Epsilon-dominance produces less uniform coverage (higher NN_CV, higher L2*) than purpose-built space-filling methods (LHS, Sobol) in the full 3D DV space. However, the constraint that solutions lie on the simplex is external to the coverage comparison. In this analytic problem, all points are Pareto-optimal, so the simplification to the hyperplane is an artificial constraint. In hydrology, the Pareto front will be much smaller and confined to the feasible drought region by physics, not by Borg's discretion.

### Recommendation

The 3D finding is expected and does NOT invalidate MOEA-FIND because:
1. In hydrology applications, the feasible drought region is unknown a priori and constrained by generator capabilities and physics. Epsilon-dominance naturally discovers this region; LHS/Sobol cannot guide sampling to it without post-hoc filtering.
2. The real advantage of MOEA-FIND emerges when the DV space is much larger (936 for 78-year traces) and the objective landscape has real trade-offs and constraints. In this regime, Borg's adaptive operator selection and epsilon-progress restart outperform blind space-filling methods.
3. Coverage uniformity is important but secondary to coverage of the feasible region. A well-distributed sample that excludes infeasible droughts is more valuable than a uniformly distributed sample that includes them.

### Next Steps

Verify coverage quality empirically when transitioning from analytic to hydrology (Experiment 2.1): compare Kirsch-based MOEA-FIND ensemble coverage vs. random Kirsch library + LHS subsampling in drought characteristic space. This will test whether the feasible-region advantage dominates over the uniformity disadvantage.

---
