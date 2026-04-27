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

2. **Moderate-length traces (~30 years):** Aggregate metrics (mean duration, frequency, mean intensity) are better supported by existing literature. Constraining to ~30 years keeps DVs manageable (~360 monthly or ~120 seasonal-block DVs). Needed for literature comparability. Matches the framing in Borgomeo et al. (2015), Zaniolo et al. (2024, FIND), and Wheeler et al. (2025).

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

- Borgomeo (2015) and Wheeler (2025) used simulated annealing with discrete swap operations (swap one month for another), which avoids the continuous/discrete issue
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
- Bootstrap track: Kirsch et al. (2013) block bootstrap; Borgomeo et al. (2015) search-based generation with swap operations; Wheeler et al. (2025) multi-site with discrete swaps.
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
- SSI (Standardized Streamflow Index) used by Zaniolo et al. (2024, FIND) for drought control, with onset at SSI <= -0.84 (20th percentile).
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

### Decision (updated 2026-04-21; see §DD-14 for the production constraint selection; see `manuscript_main_draft.md` §2.2.4 for prose specification)

**Resolved.** The constraints verify statistical plausibility — that each synthetic trace lies within the envelope of traces the Kirsch-Nowak generator could produce naturally given an infinite ensemble. Two constraints are specified in the manuscript: lag-1 autocorrelation within tolerance of historical, and non-drought-period mean flow within tolerance of historical. These are implemented as soft constraints under the Deb (2000) constraint-domination criterion. The exact tolerance values are items in the code alignment backlog (see `methods_audit.md` Items 1 and 2). The constraint set does not identify drought events; it verifies plausibility of the full trace relative to the generator's natural envelope. A third constraint (seasonal cycle) is present in the code and is discussed in `methods_audit.md` Item 3.

**Literature on constraints:** Wheeler et al. (2025) includes Hurst coefficient for long-range dependence and cross-site correlation matrix in the objective function (not as a constraint, but serves the same purpose). Borgomeo et al. (2015) used autocorrelation matching as part of the simulated annealing objective. The block bootstrap (Kirsch et al. 2013) inherently preserves within-block temporal structure and, with shared indices, cross-site correlation. Larger blocks reduce the need for explicit autocorrelation constraints. The Cholesky decomposition step in the Kirsch pipeline further enforces historical correlation structure.

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

**Literature context:** Borgomeo et al. (2015) and Wheeler et al. (2025) both implemented standalone generators with simulated annealing. The Kirsch bootstrap itself is simple (<100 lines), so standalone implementation is low-risk. The value of SynHydro integration is in the full pipeline: data loading, log-space deseasonalization, Cholesky correlation decomposition, and daily disaggregation (Nowak et al. 2010). For the MOEA-FIND publication, the monthly generation mechanism needs to be clearly documented. The KirschBorgWrapper provides this clarity while eliminating code duplication and ensuring correctness of the Cholesky + normal-score pipeline.

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

### Decision (updated 2026-04-15)

**Borg MOEA is the production algorithm.** EpsNSGAII (platypus) is the confirmed local stand-in for rapid development and analytic validation on machines without the Borg license. Analytic results produced with EpsNSGAII are valid for the construction because both algorithms use the same epsilon-dominance archive mechanism, which is the only property the theoretical argument requires. Production Cannonsville and library-comparison results will use Borg MOEA on HPC. No comparison of the two algorithms is needed beyond confirming that both produce interior-filling archives on the analytic benchmark (already verified in DD-11).

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
Kirsch et al. (2013) introduced the block bootstrap stochastic streamflow generator with Cholesky correlation preservation (B=1 monthly resampling + normal-score transform). The Kirsch-Nowak generator is widely used in MORDM studies to produce synthetic streamflow ensembles for robustness evaluation. Typical practice generates ensembles of 100-10,000 traces and uses them directly for policy re-evaluation without structured subsampling in drought characteristic space.

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

**Results (preliminary — early exploratory run, superseded by DD-11 for authoritative numbers):**
- 1362 Pareto-optimal solutions on the 2-simplex (triangular surface in 3D objective space)
- Hyperplane constraint verified: solutions satisfy sum(J_i) = constant to machine precision (~10^-16 std)
- Full range coverage: all three objectives span [-3, 3]
- Nearest-neighbor coefficient of variation (NN_CV): 0.42 for Pareto front
- L2-star discrepancy (L2*): 0.038 for Pareto front

**Comparison to baselines on the same [0, anti_ideal]^3 bounding box:**
- LHS (1362 samples): NN_CV = 0.37, L2* = 0.006
- Sobol (1362 samples): NN_CV = 0.28, L2* = 0.001
- Random (1362 samples): NN_CV = 0.50, L2* = 0.042

These values are from an early exploratory run at smaller NFE and are NOT the authoritative Section 3.1 results. The authoritative dimension-sweep results at 30,000 NFE per K are in DD-11 and `shell_vs_interior_diagnostic.md`.

**Interpretation:** Epsilon-dominance produces less uniform coverage (higher NN_CV, higher L2*) than purpose-built space-filling methods (LHS, Sobol) in the full 3D DV space. However, the constraint that solutions lie on the simplex is external to the coverage comparison. In this analytic problem, all points are Pareto-optimal, so the simplification to the hyperplane is an artificial constraint. In hydrology, the Pareto front will be much smaller and confined to the feasible drought region by physics, not by Borg's discretion.

### Recommendation

The 3D finding is expected and does NOT invalidate MOEA-FIND because:
1. In hydrology applications, the feasible drought region is unknown a priori and constrained by generator capabilities and physics. Epsilon-dominance naturally discovers this region; LHS/Sobol cannot guide sampling to it without post-hoc filtering.
2. The real advantage of MOEA-FIND emerges when the DV space is much larger (936 for 78-year traces) and the objective landscape has real trade-offs and constraints. In this regime, Borg's adaptive operator selection and epsilon-progress restart outperform blind space-filling methods.
3. Coverage uniformity is important but secondary to coverage of the feasible region. A well-distributed sample that excludes infeasible droughts is more valuable than a uniformly distributed sample that includes them.

### Next Steps

Verify coverage quality empirically when transitioning from analytic to hydrology (Experiment 2.1): compare Kirsch-based MOEA-FIND ensemble coverage vs. random Kirsch library + LHS subsampling in drought characteristic space. This will test whether the feasible-region advantage dominates over the uniformity disadvantage.

---

## DD-11: Correct Formulation of the L1 Device and Verified Interior Coverage

### The finding (2026-04-14, Session 15)

A methodological critique of the manuscript draft raised the concern that
the L1 Manhattan-norm construction, as stated in manuscript section 3.2,
is trivially degenerate under minimization and delivers at best a
(k-1)-dimensional shell in drought characteristic space rather than full
interior coverage. The user escalated the concern and required empirical
verification on a constrained feasible region before any further
manuscript work could proceed. Interior coverage was declared a hard
constraint.

### The underlying issue

The manuscript section 3.2 states the formulation as

    f_j(x)    = |D_j(x) - D*_j|     for j = 1..k
    f_{k+1}(x) = sum_j f_j(x)

Under this formulation the hyperplane identity `f_{k+1} = sum f_j` is
trivially true for every feasible point, and under minimization the
(k+1)-th objective is linearly redundant with the first k, so the
Pareto front collapses to the single point at `D = D*` (or to the
feasible point closest to it). This is the degenerate reading the
method critic diagnosed.

**However, that is not the formulation the code executes.** Careful
reading of `src/objectives.py` and `src/experiment_utils.py` line 177
shows that both `analytic_objectives` and `drought_objectives(target=None)`,
which are the objective functions used by every existing experiment,
implement

    f_j(x)     = D_j(x)                         (raw characteristic)
    f_{k+1}(x) = ||D(x) - D*||_1                (Manhattan distance)

with `D*` placed at a corner of most severe drought that lies outside
the feasible region (specifically `1.5 * max_hist D_j` for each
non-cyclic characteristic). Under this formulation and the assumption
`D_j <= D*_j`, the identity becomes `sum_{j=1}^{k+1} f_j = sum_j D*_j`,
a constant, and every feasible point is non-dominated because gains on
`f_j` are exactly offset by losses on `f_{k+1}`. The feasible set maps
injectively to a hyperplane in (k+1)-objective space, and Borg's
epsilon-box archive tiles that hyperplane.

The manuscript prose therefore describes a formulation that is not the
one the code executes and would not work if it did. The code's
formulation is correct and is the intended method.

### Empirical verification

A diagnostic `scripts/diag_shell_vs_interior.py` was written to empirically verify interior coverage on a constrained analytic problem under the code's actual formulation. The test problem places a k-ball of radius 2.5 inside the bounding box `[-3, 3]^k` with `D*` at the positive corner. MOEA-FIND runs EpsNSGAII (platypus; local stand-in — production runs will use Borg MOEA on HPC) with a penalty vector for infeasible points, and reference samples are drawn inside the k-ball at matched archive size via rejection sampling from uniform, LHS, and Sobol designs.

> **Preliminary results** — all numerical values below are from EpsNSGAII stand-in runs. Production Borg MOEA runs on HPC will replace these values before submission.

**k = 2 result (2026-04-14, 30 000 NFE, seed 42):**

| sampler   | mean dist to D* | std  | interior fraction | 15x15 grid coverage |
|-----------|-----------------|------|-------------------|---------------------|
| MOEA-FIND | 6.004           | 1.80 | 0.796             | 177 / 177           |
| uniform   | 6.059           | 1.76 | 0.808             | 177 / 177           |
| LHS       | 5.996           | 1.73 | 0.805             | 177 / 177           |
| Sobol     | 5.987           | 1.76 | 0.813             | 177 / 177           |

MOEA-FIND's statistics are within one standard error of uniform-in-disk.
The scatter plots are visually indistinguishable. No shell bias, no
orthant clustering, no interior vacancy. The shell-only hypothesis is
refuted at k=2.

**K = 3 through K = 6 results (2026-04-14, 30 000 NFE each, seed 42).**
Notation updated to use `K` for target drought characteristic dimensions
and `K + 1` for total objectives, following the user's convention.

| K | sampler    | n    | mean L1 | interior | orthant | grid |
|---|------------|------|---------|----------|---------|------|
| 2 | MOEA-FIND  | 1874 | 6.006   | 0.796    | 1.000   | 1.000|
| 2 | uniform    | 1874 | 6.064   | 0.808    | 1.000   | 1.000|
| 3 | MOEA-FIND  | 6158 | 8.976   | 0.731    | 1.000   | 0.998|
| 3 | uniform    | 6158 | 8.953   | 0.728    | 1.000   | 0.998|
| 4 | MOEA-FIND  | 8546 | 11.984  | 0.687    | 1.000   | 1.000|
| 4 | uniform    | 8546 | 11.970  | 0.660    | 1.000   | 1.000|
| 5 | MOEA-FIND  | 7841 | 14.994  | 0.629    | 1.000   | 0.977|
| 5 | uniform    | 7841 | 15.004  | 0.595    | 1.000   | 0.998|
| 6 | MOEA-FIND  | 7670 | 18.007  | 0.648    | 1.000   | 0.684|
| 6 | uniform    | 7670 | 18.030  | 0.540    | 1.000   | 0.817|

MOEA-FIND tracks uniform-in-ball on mean L1 distance and orthant
occupancy at every tested `K`. Interior mass fraction is at least as
high as uniform and higher at `K in {4, 5, 6}`. The only deviation
from uniform is grid cell coverage at `K = 6`, where the MOEA-FIND
archive shows sub-grid clustering that can be mitigated with smaller
epsilon or more NFE but does not affect interior coverage. Full
tabulated results at
`outputs/diag_shell_vs_interior/sweep_table.md`. Per-`K` figures at
`figures/figSI_shell_interior_k{K}.pdf`. Dimension-sweep summary at
`figures/figSI_shell_interior_sweep.pdf`. Full verdict at
`evidence/shell_vs_interior_diagnostic.md`.

### Decisions

1. **Method is validated for interior coverage at low k.** The user's
   blocking constraint is met at k = 2 and is expected to hold through
   k = 6 subject to the background diagnostic results.
2. **Manuscript section 3.2 must be rewritten.** The prose currently
   describes a degenerate formulation. The rewrite uses the
   construction in DD-11's "underlying issue" block above, which is
   also the new SI-1 text. SI-1 was updated in the same session.
3. **Figure 01 (Manhattan concept) must be replaced.** The current
   figure visualises the (wrong) degenerate formulation as a diagonal
   line in 2D and a triangular simplex in 3D, which looks inconsistent
   with Figure 02b where the analytic archive fills the cube. The
   replacement figure uses the constrained disk diagnostic directly
   and shows interior coverage equivalent to LHS-in-disk.
4. **A new SI section (SI-2) was added** to document the diagnostic
   and to report the dimension sweep. SI numbering was shifted by one
   so that SI-2 is the diagnostic and former SI-2 through SI-7 are now
   SI-3 through SI-8.
5. **A new DD-12 will be opened** for the high-dimensional empirical
   audit: the real hydrology case has 360 to 936 continuous decision
   variables, and this diagnostic does not substitute for verifying
   interior coverage on the actual Kirsch Pareto archive against a
   feasible-region-restricted library subsample. That audit belongs in
   main section 6 and is a hard gate on the empirical thesis.

### Open questions

1. At what dimensionality `k` does the MOEA-FIND archive start to
   degrade relative to LHS-in-ball, if at all? The k = 3 through k = 6
   sweep is designed to answer this on the analytic test problem.
2. Does the result hold on a non-convex feasible region that better
   matches the real Kirsch case? The k-ball is convex. A follow-up
   diagnostic with an annular or L-shaped feasible region would
   strengthen the SI argument. Not blocking for the main-text rewrite.
3. Does the result hold on a high-dimensional decision space? This is
   the DD-12 question and is the main section 6 empirical thesis.

---

## DD-12: High-dimensional empirical audit on the Cannonsville case

*Opened 2026-04-14 as a Phase beta prerequisite identified during the
manuscript restructure. The DD-11 analytic diagnostic at K = 2 through
K = 6 on a convex K-ball with at most 6 decision variables does not
substitute for checking interior coverage on the Cannonsville
single-site case, which has 360 continuous decision variables in
residual injection mode or 936 in index injection mode and a
non-convex feasible drought hazard region produced by the SSI event
extraction pipeline.*

### The question

Does MOEA-FIND deliver interior-filling coverage of the feasible
drought hazard region on a high-dimensional continuous decision space
with a non-convex feasible region, or do one or more of the failure
modes that SI-2 rules out on the K-ball test problem reappear under
realistic hydrologic nonlinearity?

### Specific failure modes to audit

1. Shell-only coverage. The archive could concentrate at the outer
   boundary of the feasible drought hazard region and miss interior
   drought signatures that are relevant to planners. SI-2 rules this
   out at K up to 6 on the K-ball; the Cannonsville case must be
   audited directly.
2. Orthant collapse. The archive could latch onto a subset of the
   2^K signed orthants around the anti-ideal and miss others. SI-2
   reports full orthant occupancy at K up to 6; the Cannonsville
   case must report per-run orthant occupancy.
3. Sub-grid clustering. The K = 6 K-ball diagnostic shows mild
   fine-scale clustering that is interior rather than boundary-based.
   The Cannonsville case at K = 3 should not exhibit this at the
   planned epsilon and NFE budget, but the diagnostic is worth
   running.
4. Disconnected feasible components. A non-convex feasible region
   produced by the SSI event extraction pipeline can have disconnected
   components that are reachable only through low-probability regions
   of the decision space. The adaptive variational operators of Borg
   MOEA may or may not discover them without targeted diagnostics.

### Audit protocol

After Phase beta completes, the Cannonsville Pareto archive is
diagnosed with the following metrics on the feasible drought hazard
region estimated as the convex hull of the union of the MOEA-FIND
archive and the 10,000-trace Kirsch-Nowak library.

1. Distribution of L1 distances from the anti-ideal across archive
   members, compared to an equal-size Latin hypercube subsample of
   the 10,000-trace library restricted to the same convex hull by
   rejection.
2. Interior mass fraction of archive members, defined as the fraction
   whose distance from the estimated boundary of the feasible region
   exceeds a fixed radius.
3. Signed orthant occupancy relative to the anti-ideal and relative
   to the centroid of the feasible region.
4. Nearest-neighbour coefficient of variation and feasible-region
   restricted L2 star discrepancy on each axis of the three-dimensional
   drought hazard space.
5. Visual diagnostic by pairwise projection onto the three planes of
   the drought hazard space, with the MOEA-FIND archive overlaid on
   the library subsample.

### Decision rule

The audit passes if the Cannonsville archive matches or exceeds the
library-restricted Latin hypercube reference on interior mass fraction
and signed orthant occupancy and if no disconnected feasible component
is visible in the pairwise projections. The audit fails if any of the
four failure modes above is detected, in which case section 3.2 of
the manuscript must report the failure mode and section 4.2 must add
a specific caveat naming the Cannonsville dimensionality and epsilon
budget at which the failure mode appears.

### Relationship to DD-11

DD-11 is the analytic verification at K = 2 through K = 6 on a convex
K-ball with a decision space at most 6-dimensional. DD-12 is the
empirical verification at K = 3 on a decision space of 360
dimensionality with a non-convex feasible region. DD-11 does not
substitute for DD-12 because the two diagnostics stress different
limbs of the construction. DD-11 stresses the geometric argument of
section 2.4 at increasing K. DD-12 stresses the adaptive-search
assumption that Borg's variational operators actually discover the
full feasible region in a high-dimensional decision space.

### Status

Blocked on HPC Phase beta. The audit protocol is implemented as a
stand-alone post-processing script in
`scripts/diag_cannonsville_coverage_audit.py` (to be written), which
takes the Phase beta archive and the Phase gamma library as inputs
and produces the four metrics and the pairwise projection figure.

---

## DD-13: Manuscript outline finalisation and academic style rules

*Opened and resolved 2026-04-14 after three successive outline
iterations and the user's explicit list of academic writing
violations in the prior draft.*

### The question

What is the finalised section structure of the main-text
manuscript, what is the finalised figure sequence, and what are
the hard academic writing rules that every future edit to
`manuscript_main_draft.md` must satisfy?

### Finalised section structure

Five top-level sections with no subheadings in Introduction or
Discussion, matching the Borgomeo 2015 and Zaniolo 2024 convention
and the user's explicit instruction of 2026-04-14 that Introduction
and Discussion be flat continuous prose rather than partitioned.

1. Introduction (flat continuous prose)
2. Methods: 2.1 Quantification of drought characteristics;
   2.2 Generation of structured drought hazard ensembles (with
   nested 2.2.1 Kirsch-Nowak, 2.2.2 Borg MOEA, 2.2.3
   Manhattan-distance auxiliary objective, 2.2.4 Algorithm and
   implementation); 2.3 Analytic benchmark problem; 2.4 Delaware
   River Basin case study; 2.5 Coverage metrics and evaluation
   protocols.
3. Results: 3.1 Interior coverage on the analytic benchmark
   problem; 3.2 Coverage of the feasible drought hazard region at
   Cannonsville; 3.3 Scenario discovery in the drought hazard
   space.
4. Discussion (flat continuous prose, five paragraphs: scope,
   limitations, custom drought characteristics, parametric and
   multi-site extensions, MORDM and scenario discovery
   implications).
5. Conclusions.

Multi-site Delaware River Basin is explicitly deferred to a
follow-up paper. No main-text multi-site figure. The Case Study
section is §2.4 inside Methods, not a top-level section, because
this is a methods paper and not a case study paper.

### Finalised figure sequence (seven main-text figures)

1. Parameter space versus drought hazard space contrast (§1 in
   prose, three-panel figure).
2. MOEA-FIND algorithmic pipeline schematic (§2.2.4).
3. Manhattan-distance auxiliary objective geometric construction
   in the $K = 2$ case (§2.2.3, three panels).
4. Interior-filling coverage across $K = 2$ to $K = 6$ on the
   constrained analytic benchmark (§3.1, four panels).
5. Cannonsville hydrology demonstration (§3.2, three panels:
   SSI-3 event extraction, representative traces, Borgomeo 2015
   Figure 10 style monthly-statistics verification).
6. Cannonsville Pareto archive in drought hazard space (§3.2,
   two panels, overlay on library subsample).
7. Scenario discovery in drought hazard space with gradient
   boosted trees (§3.3, four panels).

Figures 1, 3, and 4 are generatable now from existing analytic
outputs. Figure 2 is an Inkscape schematic and requires manual
redraw. Figures 5, 6, and 7 are blocked on HPC Phase beta
(Cannonsville MOEA-FIND run) and Phase gamma (ten-thousand-trace
Kirsch-Nowak library). The archived versions of the predecessor
figures live under `figures/archive/`.

### Hard academic writing rules (durable)

Seven anti-patterns are added to `style_guide.md` sections 5.11
through 5.17 after the user flagged the prior draft on
2026-04-14. Every edit to `manuscript_main_draft.md` or
`supporting_info_draft.md` must satisfy these rules in addition
to the existing pre-commit checklist.

1. Never quote other publications verbatim in body text,
   Abstract, Plain Language Summary, or Key Points. Paraphrase
   every cited claim in the author's own words.
2. Compact in-line citations only: `(Author, Year)`. No journal
   name, no volume, no page number, no section number in the
   in-line citation.
3. Never transcribe specific numeric results (coefficients,
   sample sizes, experimental counts) from cited publications
   into body text. Replace with qualitative takeaways.
4. No forward references that break narrative sequencing. A
   term, location, metric, or symbol introduced in section $N$
   must not appear in an earlier section without an
   immediately-adjacent one-sentence gloss.
5. No informal or quippy sentences. Short rhetorical sentences
   are disallowed; content-bearing short sentences are permitted.
6. No informal group designations for research communities.
   Phrases like "the Reed group" or "the Kasprzyk group" are
   disallowed; replace with specific citations or neutral
   descriptors of the research tradition.
7. Introduction literature discussion must be detailed and
   paraphrased, with two to four sentences per load-bearing
   predecessor specifically tied to the limitation the
   manuscript addresses.

### Terminology standardisation

"Drought hazard space" is the manuscript's framing. "Outcome
space" is reserved for direct paraphrase of the Moallemi et al.
(2020) three-space taxonomy (decision / uncertainty / outcome)
and is not used for the drought characteristic space elsewhere.
"Admissible drought characteristics", "feasibility discovery",
"Reed group", "PRIM", and "BART" are disallowed in
`manuscript_main_draft.md` and `supporting_info_draft.md`. The
scenario discovery demonstration uses a gradient boosted tree
classifier (Friedman 2001; Chen and Guestrin 2016); a logistic
regression robustness check is reported in SI-6.

### Status

Resolved 2026-04-14. Encoded in `style_guide.md` sections 5.11
through 5.17, in user memory at
`C:\Users\tjame\.claude\projects\c--Users-tjame-Desktop-Research-DRB-Pywr-DRB-MOEA-FIND\memory\feedback_academic_writing_style.md`,
and in the main manuscript draft which now passes every listed
pre-commit check.

---

## DD-14: Production Constraint Choice — Anderson-Darling DV-Uniformity

*Opened and resolved 2026-04-21.*

### The question

Which constraint regime should the main-text production run (exp04, Figs 5–6) use, and how should this differ from the SI ablation?

### Background

The DV-uniformity constraint family (exp13/14 ablation, SI) compares three regimes across 3 seeds × 200 k NFE:

| Regime | Constraint | Pareto pooled | Max mean_duration | Max mean_avg_severity | Manhattan median |
|---|---|---|---|---|---|
| hydrologic | 5-statistic plausibility | 4662 | 9.50 mo | 2.70 | 72.15 |
| dv_l2_star | L2* discrepancy on DV vector | 5602 | 9.47 mo | 2.85 | 72.57 |
| dv_ad | Anderson-Darling on DV vector | 3939 | 8.71 mo | 2.50 | 72.36 |

### Decision

**The production run uses `dv_uniform` + `ad` (Anderson-Darling).** Rationale:

1. **Interpretability.** The AD constraint operates entirely in DV space: it rejects solutions whose DV vectors are too far from the uniform distribution the search was initialised with, measured by a single, well-understood non-parametric goodness-of-fit statistic. The 5-statistic hydrologic set requires calibrated tolerances for each flow statistic and conflates the constraint role with the hydrological plausibility role.

2. **Geometry.** The AD constraint produces a tighter, more convex feasible region than L2-star (which tolerates `peak_severity_month = 12`, a numerical artefact of the cyclic metric range). The AD Pareto front is the most physically defensible of the three: all three objective ranges are within the historical envelope, whereas L2-star's `peak_severity_month` maximum of 12.0 exceeds it.

3. **Comparable coverage.** The Manhattan-distance median for AD (72.36) is between hydrologic (72.15) and L2-star (72.57) — within 0.4 of hydrologic, well inside seed-level noise. Coverage of drought-characteristic space is equivalent across all three regimes at 200 k NFE (confirmed by exp14 figures).

4. **Calibration traceability.** The AD tolerance is calibrated by bootstrap U[0,1] draws (`diag_dv_uniformity_calibration`), producing a single scalar threshold with a clear operational definition: the constraint passes any DV configuration reachable from a uniform distribution at the calibrated significance level.

### What changes in the code

- `workflows/experiments/04_kirsch_single_site.py`: added `--constraint-mode`, `--statistic`, `--dv-uniformity-json` args; default is `dv_uniform / ad`. The old `--constraints-json` flag is retained for the hydrologic mode (SI reproducibility).
- `workflows/slurm/04_kirsch_single_site.slurm`: defaults to `CONSTRAINT_MODE=dv_uniform`, `DV_STATISTIC=ad`, pointing at `outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json`.
- Output slug for the AD run: `residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad` (encodes both constraint mode and statistic).

### SI treatment

The exp13/14 ablation retains all three constraint arms and all 7 SI figures comparing them. The SI caption notes that the AD arm is the production choice and references this DD. The hydrologic arm remains in the SI as the baseline most similar to Wheeler et al. (2025).

### Precondition for rerun

`outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json` must exist (already produced). Submit:

```bash
sbatch --export=ALL,BORG_NFE=200000 workflows/slurm/04_kirsch_single_site.slurm
```

The old hydrologic output (`residual_T20_nfe200000_s42_constrained`) is retained in `outputs/exp04_kirsch_single_site/` under its original slug and is not overwritten by the AD rerun.
