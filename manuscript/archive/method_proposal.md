> **ARCHIVED 2026-04-15 — historical provenance only. Do not act on this document.**
> Origin design proposal. Retained for: §3.1 architecture diagram, §8 parametric reference list. Sync note at top of file identifies superseded sections. Current method state: `notes/design_decisions.md`. Current manuscript: `manuscript_main_draft.md`.

# MOEA-FIND: Method Proposal and Design Notes

*Date: 2026-04-12*
*Status: Active development*
*Origin: Trevor Amestoy's MOEA-FIND proposal presentation (April 2026)*

> **2026-04-14 sync update.** This proposal document captures the
> original design thinking and is preserved for provenance. The
> authoritative current state of the method, the manuscript
> structure, and the experiment plan is `manuscript_main_draft.md`
> together with the plan file
> `C:\Users\tjame\.claude\plans\staged-painting-dongarra.md`.
> The following points in this document have been superseded or
> refined:
>
> - Section 2, "The Manhattan Norm Trick": the formulation stated
>   here is correct in spirit but was written before DD-11. The
>   correct formulation uses $f_j = D_j$ for $j = 1, \dots, K$ and
>   $f_{K+1} = \sum_j \lvert D_j - D^*_j \rvert$ with every
>   objective minimised, and the non-dominance argument holds on
>   the image of the generator map rather than in the decision
>   space. See `manuscript_main_draft.md` Section 2.2.3 and SI-1
>   for the finalised exposition.
> - Section 2.2, "Why This Works": the shell-only concern has been
>   empirically refuted at $K = 2$ through $K = 6$ on a constrained
>   $K$-ball benchmark. See DD-11 and
>   `shell_vs_interior_diagnostic.md`.
> - Section 3.2, "Multi-Site Generation": the multi-site Delaware
>   River Basin extension is no longer in the main-text scope of
>   the present paper and is deferred to follow-up. The manuscript
>   demonstrates on the single-site Cannonsville inflow only. Plan
>   file Phase 4 and Phase 5.6 reflect this change.
> - Section 5, "Research Questions": see `research_questions.md`
>   for the updated and synced version. The 2026-04-14 sync header
>   in that file documents the specific corrections (Zaniolo 2024,
>   Wheeler 2025, no BART or PRIM, multi-site deferred, hazard
>   space framing).
> - Section 7, "Experimental Plan": Phase 3 parametric CDF
>   formulation and Phase 4 multi-site DRB are both deferred to
>   follow-up papers. Phase 1 analytic and Phase 2 Kirsch
>   single-site are the scope of the present paper, corresponding
>   to main-text Sections 3.1 and 3.2 through 3.3 respectively.
> - Terminology: "outcome space", "feasibility discovery",
>   "admissible drought characteristics", "Reed group", and
>   "BART scenario discovery" are disallowed in any manuscript
>   prose per style guide sections 5.11 through 5.17. The
>   downstream scenario discovery demonstration in main-text
>   Section 3.3 uses a gradient boosted tree classifier. Any
>   text in this proposal file that contains those terms is
>   historical and should not be transcribed into the manuscript.

---

## 1. Problem Statement

Robust evaluation of water supply operating policies requires testing under diverse drought scenarios that span the space of plausible drought characteristics (frequency, intensity, duration). Existing approaches have two limitations:

1. **Library-and-subsample:** Generate a large library of synthetic traces (e.g., 10,000 via Kirsch-Nowak), characterize each by drought metrics, then subsample using LHS or KD-tree to get representative coverage. This is indirect: coverage quality depends on the library being large enough and the subsampling algorithm.

2. **Single-objective search (FIND, Borgomeo):** Use optimization to find a synthetic trace matching one target drought characteristic. Generates one scenario per optimization run. To build an ensemble, run many times with different targets, with no guarantee of structured coverage.

MOEA-FIND solves both problems by formulating scenario generation as a multi-objective optimization problem where the Pareto front IS the structured ensemble.

---

## 2. Core Innovation: The Manhattan Norm Trick

### 2.1 Setup

Given k drought characteristics we want to vary (e.g., frequency, duration, intensity), define k+1 objectives:

- J_1 = D_1 (e.g., relative drought frequency)
- J_2 = D_2 (e.g., relative drought duration)
- ...
- J_k = D_k
- J_{k+1} = ||D - D*||_1 (L1 Manhattan norm from anti-ideal point D*)

where D* is placed at the corner of objective space opposite the ideal point (i.e., at the maximum plausible values of all drought characteristics).

### 2.2 Why This Works

**Claim:** The Pareto front of this (k+1)-objective problem lies on a hyperplane in objective space.

**Intuition (2D case):** With J_1 = X_1, J_2 = X_2, J_3 = |X_1 - X*_1| + |X_2 - X*_2|:
- When minimizing all three, for any point (X_1, X_2), J_3 is determined by J_1 and J_2
- If X* is at (max, max), then J_3 = (X*_1 - X_1) + (X*_2 - X_2) for all Pareto-optimal points
- This means J_1 + J_2 + J_3 = X*_1 + X*_2 = constant
- The Pareto front lies on the plane J_1 + J_2 + J_3 = C (a simplex)

**Generalization to k dimensions:** For k drought objectives plus Manhattan norm, the Pareto front lies on the (k-1)-simplex defined by sum(J_i) + J_{k+1} = sum(D*_i). Borg's epsilon-dominance archiving tiles this simplex with approximately uniform epsilon-box grid spacing.

**Result:** The Pareto front, projected back to the k-dimensional drought characteristic space, provides near-uniform coverage. This was demonstrated in the 2D proof-of-concept (slides 13-14) where Borg's Pareto front was visually indistinguishable from an LHS sample.

### 2.3 Proof-of-Concept Implementation (2D)

From the proposal presentation (slides 7-15):

**Setup:**
- Decision variables: X_1, X_2 (continuous, range [-3, 3])
- Objectives: J_1 = X_1, J_2 = X_2, J_3 = ||X - X*||_1
- Ideal point: (-3, -3) (minimum of both)
- Anti-ideal point X*: (3, 3) (maximum of both)
- J_3 = |X_1 - 3| + |X_2 - 3| = (3 - X_1) + (3 - X_2) = 6 - X_1 - X_2

**Results:**
- Sobol samples lie on the Pareto front (the 3D surface is a tilted plane)
- Borg generates ~8,365 Pareto-optimal solutions
- Projected to (X_1, X_2), the Borg solutions are nearly indistinguishable from a Latin hypercube sample
- Both checklist items verified: (1) structured samples lie on PF, (2) PF projection reflects structured samples

---

## 3. Practical Application: Synthetic Drought Generation

### 3.0 Two Sampling Formulations

The method can use two different mechanisms for generating candidate monthly flows:

**Formulation A: Block Bootstrap (Kirsch-Nowak)**
Decision variables are bootstrap sample indices. Generated flows are recombinations of historical observations. Temporal/spatial correlation is preserved by shared block indices. Limited to the historical envelope of observed flows.

**Formulation B: Parametric CDF + Vine Copula**
Decision variables are CDF probabilities p in [0,1]. Monthly flows are generated by inverse CDF from fitted parametric distributions (kappa distribution recommended, Svensson et al. 2017). Temporal dependence imposed via periodic D-vine copula (lag-1, lag-2). Spatial dependence via multivariate vine. This allows generating flows outside the historical range, enabling exploration of more extreme drought characteristics.

Both formulations use the same Borg MOEA + Manhattan norm framework. The comparison between them (does parametric sampling meaningfully expand the achievable drought space?) is itself a contribution.

### 3.1 Architecture (slides 16-25)

```
Q_H (historical) ──→ Calc Reference Drought Stats ──→ D_ref
                                                          │
                                                          ▼
Borg MOEA ──→ vars (DV indices) ──→ KirschBorgWrapper ──→ Q_Si
    ▲                                 (SynHydro pipeline)   │
    │                                                       ▼
    └──── J (objectives) ◄──── Calc Drought Stats ◄──── D_Si
           C (constraints)
```

**Decision variables:** Array of monthly or yearly sample indices, size varies by trace length and DV formulation.
Two formulations supported by the KirschBorgWrapper:
1. **Index injection:** Each DV is a continuous [0,1] value mapped to a historical year index. The full Kirsch pipeline (normal-score transform, Cholesky decomposition, Dec-Jan handling) ensures temporal correlation preservation.
2. **Residual injection:** Each DV is a CDF probability mapped to a standardized residual. Residuals are passed through SynHydro's Cholesky + normal-score pipeline to impose historical correlation structure.

**Objectives (for 3-characteristic case):**
- J_1 = D_R,Freq = (drought frequency of Q_Si) / (drought frequency of Q_H)
- J_2 = D_R,Dur = (mean drought duration of Q_Si) / (mean drought duration of Q_H)
- J_3 = D_R,Sev = (mean average severity of Q_Si) / (mean average severity of Q_H)
- J_4 = ||D_R - D*_R||_1 (Manhattan norm from anti-ideal)

**Constraints (slide 22):**
- Autocorrelation of generated trace must be within acceptable range of historical
- Non-drought period statistics (mean, variance) must be plausible
- Purpose: prevent Borg from finding degenerate traces (e.g., all-drought) that hit target metrics

### 3.2 Multi-Site Generation

The Kirsch block bootstrap naturally handles multi-site correlation by using the same sampling indices across all sites (slide 26). This is a key advantage over FIND, which was developed for single-site generation.

For the DRB case study, sites include:
- Cannonsville inflow
- Pepacton inflow
- Neversink inflow
- Lateral inflows along the Delaware mainstem

### 3.3 Demand Time Series

Three options (slide 26):
1. **Treat demand as another "site"** in the bootstrap. Correlation imposed automatically via shared indices. Simplest approach.
2. **Explicit stochastic model** with specified correlation to streamflow. More flexible but more complex.
3. **Regression:** f(streamflow anomalies) -> demand anomalies (done in Hadjimichael et al.). Preserves physical relationship.

Option 1 is recommended for simplicity and consistency.

---

## 4. Generalization to Higher Dimensions

### 4.1 Theoretical Scaling

For k drought characteristics + Manhattan norm = k+1 objectives:
- k=2 (3 obj): Pareto surface is a line segment in 2D drought space. Proven in POC.
- k=3 (4 obj): Pareto surface is a triangle (2-simplex) in 3D drought space. Borg handles 4 objectives routinely.
- k=4 (5 obj): Pareto surface is a tetrahedron (3-simplex). 5 objectives is within Borg's many-objective sweet spot.
- k=5+ (6+ obj): Epsilon-box coverage becomes sparser. Diminishing returns on uniformity.

**Practical recommendation:** k=3 (frequency, duration, intensity) is the natural choice for the DRB study. This gives 4 total objectives, well within Borg's capabilities.

### 4.2 Open Questions for Generalization

1. **Does the uniformity property hold in 3+ drought dimensions?** The 2D POC is compelling but needs verification in higher dimensions. Planned experiment: 3D analytic test (X1, X2, X3, Manhattan norm).

2. **Epsilon values:** How to set epsilon for each drought objective? In 2D, a single epsilon produces ~8,365 solutions. In 3D with the same epsilon, the simplex has more volume to tile but each epsilon box is also 3D. The number of Pareto solutions scales as O(1/epsilon^{k-1}).

3. **Anti-ideal placement:** In 2D, D* at (max, max) works. In higher dimensions, the anti-ideal determines the shape of the sampled region. Should it be placed at the maximum observed historical values? At some multiple thereof?

4. **Coverage verification:** How to quantify that the Pareto front provides "near-uniform" coverage? Discrepancy measures (star discrepancy, L2 discrepancy) or comparison to LHS/Sobol baselines.

---

## 5. Research Questions

### Primary

**RQ1:** Can multi-objective evolutionary optimization with an L1 norm objective produce ensembles of synthetic streamflow traces with structured, near-uniform coverage of drought characteristic space?

**RQ2:** How does the coverage quality (uniformity, extent) scale with the number of drought characteristics (dimensionality of the objective space)?

**RQ3:** Does coupling Borg MOEA with the Kirsch-Nowak block bootstrap generator produce physically plausible traces while achieving target drought characteristics?

### Secondary

**RQ4:** How does MOEA-FIND compare to existing approaches (library subsampling, FIND, Borgomeo single-objective) in terms of coverage quality, computational efficiency, and trace plausibility?

**RQ5:** What is the effect of epsilon values, anti-ideal placement, and constraint formulation on ensemble quality?

**RQ6:** Can the method be applied to multi-site basins with preserved spatial correlation?

---

## 6. Outstanding Issues and Design Decisions

### 6.1 Decision Variable Formulation

**Discrete indices vs. CDF probabilities:**
- Discrete: each DV is an integer in [1, N_years_historical]. Simple conceptually but creates a discrete search space that Borg's operators (SBX, DE, etc.) may struggle with.
- CDF: each DV is p in [0,1], mapped to a historical month via inverse CDF. Continuous search space, better for Borg. But the mapping introduces quantization artifacts.
- **Recommendation:** CDF probability formulation. Borg's operators are designed for continuous variables. The Kirsch generator can accept probabilities and map to months internally.

### 6.2 Block Size and Correlation Imposition (RESOLVED 2026-04-13)

**Correction:** The Kirsch method uses B=1 (independent monthly resampling), then imposes temporal correlation post-hoc via Cholesky decomposition of the 12x12 monthly correlation matrix in normal-score space. This is distinct from a general block bootstrap where B>1 preserves within-block autocorrelation mechanically. The two approaches are alternatives, not the same method.

**Kirsch pipeline (B=1 + Cholesky):**
1. Resample individual months independently (B=1)
2. Normal-score transform the standardized residuals
3. Apply Cholesky factor of the historical monthly correlation matrix
4. Invert the normal-score transform
5. Handle Dec-Jan transitions via a shifted "Y_prime" matrix with its own Cholesky

This produces autocorrelation structure that matches the historical record without relying on block size. Analysis revealed that earlier implementations were missing steps 2-5, which caused autocorrelation and flow duration curve anomalies in proof-of-concept experiments.

**Implementation decision:** Use the KirschBorgWrapper (see DD-06) which delegates to SynHydro's validated Kirsch implementation. This ensures all five pipeline steps are correctly executed. SynHydro's implementation is both validated and publication-ready.

**Alternative approach (not recommended for publication):** Mechanical block bootstrap with B>1 can preserve within-block temporal structure without Cholesky, but introduces block-boundary discontinuities and limits recombination flexibility. This is simpler to implement but less flexible and produces lower-quality traces. Should not be labeled as "Kirsch method."

### 6.3 Constraint Design

Constraints enforce trace plausibility (slide 22). Candidates:
1. **Lag-1 autocorrelation:** |rho_1(Q_S) - rho_1(Q_H)| < threshold
2. **Annual mean/variance:** within X% of historical for non-drought years
3. **Seasonal cycle:** monthly means within acceptable range
4. **Hurst coefficient:** long-range dependence preserved (from Wheeler et al.)
5. **Cross-site correlation:** correlation matrix of Q_S within tolerance of Q_H (for multi-site)

The block bootstrap with shared indices inherently preserves much of this structure. Constraints may be light-touch rather than aggressive.

### 6.4 Drought Metric Definitions

Need precise definitions for the objectives:
- **Threshold:** Monthly flow below P20 (20th percentile) of historical? Or a fixed absolute threshold?
- **Frequency:** Number of drought events per decade
- **Duration:** Mean number of consecutive months below threshold
- **Intensity:** Mean deficit below threshold during drought months (cumulative or average?)

These must match the definitions used in the policy re-evaluation to ensure consistency.

### 6.5 Computational Cost

Each Borg evaluation = run Kirsch generator + compute drought stats. If the generator processes a 78-year trace:
- Kirsch generator: ~0.01-0.1s (block bootstrap, fast)
- Drought stat calculation: ~0.001s (simple threshold operations)
- Total per evaluation: ~0.01-0.1s

With 936 DVs and 4 objectives, Borg may need 500K-1M NFE for good coverage. At 0.1s/eval, that's ~14-28 hours serial. Parallelizable via Borg's master-worker architecture.

### 6.6 SynHydro Integration (RESOLVED 2026-04-13)

MOEA-FIND wraps SynHydro's validated Kirsch-Nowak pipeline via the KirschBorgWrapper (see DD-06). The wrapper handles:
1. DV-to-index or DV-to-residual mapping (two injection modes)
2. Passing mapped DVs to SynHydro's Kirsch generator
3. Receiving generated synthetic traces
4. Computing drought statistics on traces
5. Interfacing with Borg's Python API

SynHydro exposure requirement: The KirschGenerator class must expose a method accepting externally provided year indices or standardized residuals, running the full Cholesky + normal-score pipeline, and returning synthetic flows. This is a small, clean addition to the existing class and ensures MOEA-FIND does not duplicate the complex Cholesky and Dec-Jan boundary handling logic.

---

## 7. Experimental Plan

### Phase 1: Proof-of-Concept Replication and Extension

**Experiment 1.1: COMPLETE (2026-04-13)**
Replicated 2D analytic test (J_1=X_1, J_2=X_2, J_3=Manhattan). Verified near-uniform coverage. Computed L2-star discrepancy and nearest-neighbor CV. Compared to LHS and Sobol baselines. *Literature basis: Hadka and Reed (2013) for Borg MOEA configuration; IEEE CEC (2019) work on sampling reference points on Pareto fronts of benchmark problems.*

**Experiment 1.2: COMPLETE (2026-04-13)**
Extended to 3D analytic test (J_1=X_1, J_2=X_2, J_3=X_3, J_4=Manhattan). Verified that 1362 Pareto solutions satisfy hyperplane constraint (sum(J_i)=C) to machine precision. Coverage metrics: NN_CV=0.42 (Pareto) vs 0.37 (LHS) vs 0.28 (Sobol) vs 0.50 (Random); L2*=0.038 (Pareto) vs 0.006 (LHS) vs 0.001 (Sobol) vs 0.042 (Random). Conclusion: epsilon-dominance produces less uniform coverage than purpose-built space-filling methods in full DV space, but this is expected because all points are Pareto-optimal in the analytic test. Real advantage emerges in constrained feasible regions (hydrology). See DD-10 for full analysis. *Literature basis: Deb and Jain (2014) NSGA-III for reference-point coverage comparison; Zhang and Li (2007) MOEA/D for decomposition-based comparison.*

**Experiment 1.3: PENDING**
Sensitivity analysis on epsilon values and NFE budget. How does coverage quality (discrepancy) depend on these parameters? Number of Pareto solutions scales as O(1/epsilon^{k-1}). *Literature basis: Hadka and Reed (2013) epsilon-progress restart triggering and adaptive operator selection.*

### Phase 2: Synthetic Streamflow Generation (Single-Site, Kirsch Bootstrap)

**Experiment 2.1:** Single-site Kirsch-based generation via KirschBorgWrapper with 2 drought objectives (duration, mean average severity) + Manhattan. Verify that generated traces have the expected drought characteristics. The Kirsch pipeline operates in log-space with Cholesky decomposition (normal-score transform, correlation matrix factor, inverse transform, Dec-Jan handling) to preserve historical autocorrelation (Kirsch et al. 2013). Drought threshold defined as P20 monthly flow following SSI convention (onset at SSI <= -0.84, 20th percentile). *Literature basis: Kirsch et al. (2013) for Kirsch-Nowak methodology; Borgomeo et al. (2015) for search-based generation precedent; Svensson et al. (2017) for drought threshold selection and SDF curve context.*

**Experiment 2.2:** Add plausibility constraints. Assess trade-off between constraint tightness and coverage quality. Candidate constraints: lag-1 autocorrelation tolerance, annual statistics within range, Hurst coefficient for long-range dependence. The Kirsch pipeline inherently preserves much temporal structure via Cholesky decomposition, so constraints may be light-touch. *Literature basis: Wheeler et al. (2025) for Hurst coefficient in search-based generation; Kirsch et al. (2013) for Cholesky pipeline; HHJ and NPPI methods for optimal block length estimation.*

**Experiment 2.3:** Compare coverage to: (a) random Kirsch library + LHS subsampling, (b) repeated single-objective FIND runs with designed target grid (Zaniolo et al. 2024), (c) Sobol sequence in drought space with nearest-neighbor trace selection, (d) random uniform Kirsch library (unstructured baseline). Coverage quality measured by L2 discrepancy, star discrepancy, NN_CV, and MST mean (following Bonham et al. 2024). This is the core comparison for RQ4. *Literature basis: Zaniolo et al. (2024) FIND for single-objective comparison; Borgomeo et al. (2015) for single-target generation; LHS theory for space-filling design comparison; Sobol sequences for quasi-random low-discrepancy baseline.*

**Experiment 2.4:** Event-level traces (1-8 year events) with event-specific objectives (event duration, peak intensity, cumulative severity). Compare achievable drought space to Experiment 2.1 (aggregate metrics over 30-year traces). This addresses DD-01 and demonstrates both temporal framings. The event-level approach is the more novel contribution (see Critical Analysis, Section 1). *Literature basis: Climate stress testing literature on splicing drought events into baseline sequences; SDF (Severity-Duration-Frequency) curves for joint drought characterization; VLB (Variable-Length Block) methods for generating flows outside historical record.*

### Phase 3: Parametric CDF Formulation

**Experiment 3.1:** Fit four-parameter kappa distributions to historical monthly flows at each site. Validate goodness-of-fit (KS test, QQ plots). Kappa distribution is bounded below at zero with flexible tail behavior for low-flow extremes. *Literature basis: Svensson et al. (2017) established kappa as the flexible choice for monthly streamflow with superior low-flow characterization.*

**Experiment 3.2:** Implement D-vine copula for temporal dependence (lag-1, lag-2). Validate autocorrelation structure of generated traces against historical. Consider periodic vine copula to account for seasonally varying dependence. *Literature basis: Brechmann et al. (2017) PAR(p)-vine copula for stochastic streamflow; Li et al. (2023) mixed D-vine conditional quantile model; periodic spatial vine copulas (2019, Water Resour. Manage.) for seasonal dependence.*

**Experiment 3.3:** Single-site parametric CDF + copula with same drought objectives as Phase 2. Compare achievable drought space to bootstrap formulation. Key question: does parametric sampling meaningfully extend the range of achievable drought characteristics? Also consider: (a) KDE-smoothed CDF as middle ground (mixed copula-KDE framework, 2025), (b) phase randomization (Papalexiou 2019) as alternative preserving spectral properties. *Literature basis: Papalexiou and Koutsoyiannis for phase randomization and CoSMoS framework; Mixed copula-KDE (2025, Stoch. Environ. Res. Risk Assess.); BEUM bootstrapped CDF uncertainty approach (2022).*

### Phase 4: Multi-Site DRB Application

**Experiment 4.1:** Multi-site (4 DRB inflow sites) with shared bootstrap indices. Verify spatial correlation preservation via cross-site correlation matrix comparison. The shared-index approach is simpler than Wheeler et al.'s explicit cross-correlation objective but should achieve similar results for sites within the same basin. *Literature basis: Wheeler et al. (2025) for multi-site generation with cross-correlation objective; Nowak et al. (2010) for multi-site disaggregation with K-NN conditioning; Kirsch et al. (2013) for shared-index spatial correlation.*

**Experiment 4.2:** 3 drought objectives (frequency, duration, intensity) + Manhattan = 4 objectives. Generate structured drought ensemble for DRB. This is the core deliverable: a structured ensemble with near-uniform coverage of 3D drought space, suitable for MORDM-style policy re-evaluation. *Literature basis: Herman et al. and Reed Group MORDM applications; OpenMORDM workflow documentation; Hadjimichael et al. (2020) for stakeholder-contingent robustness requiring diverse scenario coverage.*

**Experiment 4.3:** Re-evaluate NYCOptimization Pareto policies on the MOEA-FIND ensemble. Compare vulnerability analysis results to library-subsampling baseline. Use BART for scenario discovery to identify drought characteristics driving policy failure. Robustness metrics should be stakeholder-specific following Hadjimichael et al. (2020) and Sunkara et al. (2023). *Literature basis: EMODPS framework (Reed group) for policy-scenario pairing; Chipman et al. (2010) BART for scenario discovery; Hadjimichael et al. (2020) and Sunkara et al. (2023) for stakeholder-contingent vulnerability analysis.*

### Phase 5: Paper Figures and Analysis

- **Coverage quality comparison** (MOEA-FIND vs. LHS vs. FIND vs. random): L2 discrepancy and star discrepancy metrics, visual comparison of projected Pareto fronts. Baseline comparisons to Sobol sequences for quasi-random reference.
- **Bootstrap vs. parametric CDF:** Achievable drought space comparison, with emphasis on tail behavior and whether kappa + vine copula meaningfully extends beyond historical envelope.
- **Event-level vs. aggregate-level:** Compare 1-8 year event traces (novel contribution) to 30-year aggregate traces (literature-standard). Assess what each reveals about drought space structure and policy vulnerability.
- **Scaling analysis** (2D through 5D drought characteristics): Empirical verification of coverage quality on k-simplices for k=2,3,4,5. Addresses RQ2.
- **DRB case study results:** Structured ensemble generation, policy re-evaluation, BART scenario discovery. Demonstrates practical value for MORDM.
- **Trace plausibility assessment:** Autocorrelation, seasonal cycle, spatial correlation, Hurst coefficient. Comparison of bootstrap vs. parametric traces.

---

## 8. Key References

**Direct predecessors (search-based generation):**
- Borgomeo, E., et al. (2015). Numerical rivers: A synthetic streamflow generator for water resources vulnerability assessments. WRR, 51(7). — *First search-based streamflow generation (single-objective, simulated annealing)*
- Zaniolo, M., Fletcher, S., and Mauter, M. (2023). FIND: A synthetic weather generator to control drought Frequency, INtensity, and Duration. Environ. Model. Softw. — *Single-objective FID control via SSI, simulated annealing, MATLAB*
- Wheeler, K.G., et al. (2024). Multisite nonparametric stochastic streamflow generation for the Eastern Nile Basin. J. Hydrol. Eng., 30(1). — *Multi-site SA with cross-correlation objective and Hurst coefficient*

**MOEA methodology:**
- Hadka, D., and Reed, P. (2013). Borg: An auto-adaptive many-objective evolutionary computing framework. Evol. Comput., 21(2). — *Epsilon-dominance archiving, epsilon-progress restart, adaptive operator selection*
- Deb, K., and Jain, H. (2014). An evolutionary many-objective optimization algorithm using reference-point-based nondominated sorting. IEEE Trans. Evol. Comput. — *NSGA-III, reference-point coverage comparison*
- Zhang, Q., and Li, H. (2007). MOEA/D: A multi-objective evolutionary algorithm based on decomposition. IEEE Trans. Evol. Comput. — *Decomposition-based comparison*

**Kirsch-Nowak generator:**
- Kirsch, B.R., et al. (2013). Evaluating the impact of alternative hydro-climate scenarios on transfer agreements. J. Water Res. Plan. Manage., 139(4), 396-406. — *Block bootstrap in log-space, Cholesky decomposition, shared indices*
- Nowak, K., et al. (2010). A nonparametric stochastic approach for multisite disaggregation of annual to daily streamflow. WRR, 46, W08529. — *K-NN conditioning, multi-site disaggregation*

**Parametric generation and copulas:**
- Svensson, C., et al. (2017). Statistical distributions for monthly aggregations of precipitation and streamflow in drought indicator applications. WRR, 53(12). — *Four-parameter kappa distribution for monthly flows*
- Brechmann, E.C., et al. (2017). PAR(p)-vine copula based model for stochastic streamflow scenario generation. Stoch. Environ. Res. Risk Assess. — *Periodic autoregression + vine copula*
- Li, C., et al. (2023). Mixed D-vine copula-based conditional quantile model for stochastic monthly streamflow simulation. J. Hydro-environ. Res. — *D-vine with lag-1/lag-2*
- Papalexiou, S.M. (2019/2022). Phase randomization and CoSMoS framework. WRR/HESS. — *Spectral-preserving generation, flexible distribution modeling*

**Drought characterization:**
- McKee et al. — SPI; Vicente-Serrano et al. — SPEI; SSI (Standardized Streamflow Index) — *Drought index hierarchy*
- SDF curves (bivariate copula-based frequency analysis) — *Joint drought magnitude-duration characterization*

**Library-and-subsample approach (comparison baseline):**
- Bonham, N., Kasprzyk, J., Zagona, E., and Rajagopalan, B. (2024). Subsampling and space-filling metrics to test ensemble size for robustness analysis with a demonstration in the Colorado River Basin. Environ. Model. Softw., 172, 105933. — *cLHS subsampling from 500-scenario library in uncertainty input space; space-filling metrics (MSTmean, mindist) predict ranking accuracy (R2=0.77-0.91). Closest published precedent, but operates in input space not drought outcome space*
- McKay, M.D., Beckman, R.J., and Conover, W.J. (1979). A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. Technometrics, 21(2), 239-245. — *Original LHS reference*
- Sobol, I.M. (1967). On the distribution of points in a cube and the approximate evaluation of integrals. USSR Computational Mathematics and Mathematical Physics, 7(4), 86-112. — *Sobol quasi-random sequences*
- Pianosi, F., et al. (2016). Sensitivity analysis of environmental models: A systematic review with practical workflow. Environ. Model. Softw., 79, 214-232. — *LHS and Sobol as standard space-filling designs*
- Note: Bonham et al. (2024) subsample in uncertainty input space for computational efficiency. The specific workflow of subsampling in drought *characteristic/outcome* space for structured coverage does not appear in published literature. It is constructed here as a natural baseline combining standard DOE methods (LHS/Sobol) with standard generation (Kirsch-Nowak).

**Existing ensemble usage in MORDM (for contrast):**
- Hadjimichael, A., et al. (2020). Defining robustness, vulnerabilities, and consequential scenarios. Earth's Future, 8(6). — *LHS in generator parameter space (HMM transitions, mean/variance multipliers), not post-hoc subsampling in drought space*
- Quinn, J.D., et al. (2018). Exploring how changing monsoonal dynamics and human pressures challenge multi-reservoir management. EJOR, 270, 1191-1204. — *CMIP5 projections + bias-corrected downscaling, not synthetic library subsampling*
- Herman, J.D., et al. (2016). Synthetic drought scenario generation to support bottom-up water supply vulnerability assessments. JWRPM, 142(11). — *SSI-based drought control during generation, not post-hoc selection*

**Coverage metrics and sampling theory:**
- L2 discrepancy and star discrepancy — *Coverage quality metrics (Koksma-Hlawka inequality)*
- IEEE CEC (2019) — *Pareto front as space-filling design, SPREAD algorithm*

**Policy re-evaluation and scenario discovery (consumers of MOEA-FIND output):**
- EMODPS framework (Reed group) — *Policy optimization + scenario re-evaluation workflow*
- Chipman, H.A., et al. (2010). BART: Bayesian additive regression trees. Ann. Applied Stat., 4(1). — *Modern scenario discovery replacing PRIM*
- Hadjimichael, A., et al. (2020). Defining robustness, vulnerabilities, and consequential scenarios. Earth's Future, 8(6). — *Stakeholder-contingent robustness*
- Sunkara, et al. (2023). How should diverse stakeholder interests shape evaluations of complex water resources systems robustness. Earth's Future, 11(7). — *BART scenario discovery with diverse robustness definitions*

**Scenario discovery (downstream analysis):**
- BART for probabilistic scenario discovery (Water Programming blog, 2024)
- Chipman, H.A., et al. (2010). BART: Bayesian additive regression trees. Ann. Appl. Stat.

**Internal variability and climate uncertainty:**
- Lehner, F., and Deser, C. (2023). Origin, importance, and predictive limits of internal climate variability. — *Motivation: at regional scales, internal variability dominates (slide 2)*
