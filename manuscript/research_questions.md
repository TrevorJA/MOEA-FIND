# Research Questions

> **2026-04-14 sync update.** This note predates the finalised
> manuscript outline. The authoritative current state is captured in
> `manuscript_main_draft.md` and in the plan file at
> `C:\Users\tjame\.claude\plans\staged-painting-dongarra.md`. Key
> alignment points: (a) the manuscript now uses a five-section flat
> outline (Introduction, Methods, Results, Discussion, Conclusions)
> with no Introduction subheadings; (b) the hazard-space framing
> ("drought hazard space" and "feasible drought hazard region")
> replaces earlier "outcome space" and "feasibility discovery"
> language; (c) the multi-site Delaware River Basin extension is
> explicitly deferred to follow-up work and is not a main-text
> result (see RQ7 below, which is retained as a research direction
> rather than a paper-scope commitment); (d) the scenario discovery
> demonstration in Section 3.3 uses a gradient boosted tree
> classifier (Friedman 2001; Chen and Guestrin 2016), not BART or
> PRIM (see RQ8 below, which should be read with this correction);
> (e) the shell-versus-interior question raised in RQ2 has been
> resolved empirically at $K = 2$ through $K = 6$ and is recorded
> as DD-11 in `design_decisions.md`; (f) citations corrected:
> Zaniolo et al. (2024), Wheeler et al. (2025).

## Primary (Paper-Level)

**RQ1: Structured Coverage via Multi-Objective Search**
Can a multi-objective evolutionary algorithm (Borg MOEA) with an L1 norm auxiliary objective produce synthetic streamflow ensembles with near-uniform coverage of a user-defined drought characteristic space?

*How to evaluate:* L2 discrepancy and star discrepancy metrics compared to LHS and Sobol sequence baselines. The Pareto front, projected to drought characteristic space, should achieve discrepancy comparable to or better than designed space-filling samples. The Manhattan norm forces solutions onto a (k-1)-simplex where epsilon-dominance naturally tiles (IEEE CEC 2019 work on Pareto front sampling; SPREAD algorithm).

*Key references:* Hadka and Reed (2013) for Borg epsilon-dominance archiving; LHS theory for space-filling baseline; Sobol sequences for quasi-random low-discrepancy baseline; Koksma-Hlawka inequality for discrepancy bounds.

**RQ2: Dimensionality Scaling**
How does coverage quality (uniformity, extent, epsilon-box filling) scale as the number of drought characteristics increases from 2 to 5?

*How to evaluate:* Analytic tests in 2D (line), 3D (triangular simplex), 4D (tetrahedral simplex), 5D. Number of unique epsilon-boxes on a k-simplex scales as O(C^k / epsilon^{k-1}). In high dimensions, most boxes may be empty since the simplex has measure zero in the full objective space (Critical Issue 3).

*Key references:* Deb and Jain (2014) NSGA-III reference-point coverage as comparison (different mechanism); Zhang and Li (2007) MOEA/D decomposition as comparison.

**RQ3: Physical Plausibility**
Can the block bootstrap formulation, coupled with lightweight plausibility constraints, produce traces that are statistically indistinguishable from historical records in non-drought properties (autocorrelation, seasonal cycle, spatial correlation)?

*How to evaluate:* Compare lag-1 autocorrelation, seasonal means, Hurst coefficient, and cross-site correlation matrices between generated and historical traces. The Kirsch block bootstrap (Kirsch et al. 2013) with shared indices and Cholesky decomposition inherently preserves much of this structure.

*Key references:* Kirsch et al. (2013) for bootstrap plausibility; Wheeler et al. (2025) for Hurst coefficient and cross-correlation preservation; HHJ and NPPI methods for block size selection affecting plausibility.

## Secondary (Method Design)

**RQ4: Comparison to Alternatives**
How does MOEA-FIND compare to (a) LHS subsampling from a large Kirsch-Nowak library, (b) repeated single-objective FIND runs (Zaniolo et al. 2023), and (c) Borgomeo-style simulated annealing (Borgomeo et al. 2015), in terms of coverage quality and computational cost?

*Key advantages to demonstrate:* MOEA-FIND discovers the feasible drought region automatically (FIND fails silently at infeasible targets); produces structured coverage without manual target specification; single optimization run vs. N independent runs; Pareto front has mathematical relationship to drought space (simplex structure).

*Key references:* Zaniolo et al. (2024) FIND for single-objective FID control; Borgomeo et al. (2015) for single-target SA generation; Wheeler et al. (2025) for multi-site SA extension.

**RQ5: Short vs. Long Traces**
Is it more effective to generate short event-focused traces (1-8 years, 1-2 droughts with event-specific metrics) or long multi-event traces (~30 years with aggregate metrics)?

*How to evaluate:* Compare achievable drought space extent, DV tractability (60-96 vs. 360 monthly DVs), interpretability, and Borg convergence speed. Event-level metrics (duration, peak intensity, cumulative severity) vs. trace-level metrics (mean duration, frequency, mean intensity). SDF curves provide the framework for event-level characterization.

*Key references:* Climate stress testing literature on splicing events into baseline sequences; Borgomeo et al. (2015) and FIND both use moderate-length traces with aggregate metrics; SDF bivariate frequency analysis for event-level context.

**RQ6: Epsilon and Anti-Ideal Sensitivity**
How do epsilon values and anti-ideal point placement affect ensemble properties?

*How to evaluate:* Systematic sweep of epsilon values (controlling ensemble size and spacing) and anti-ideal placement (controlling extent of sampled drought space). If anti-ideal is too far from the actual Pareto front, epsilon-dominance becomes less effective; if too close, the front may not span enough space. Epsilon-progress restart triggering (Hadka and Reed 2013) helps escape local optima as epsilon changes.

## Applied (DRB Case Study)

**RQ7: Multi-Site Basin Application**
Can MOEA-FIND generate structured drought ensembles for the multi-site Delaware River Basin while preserving inter-site correlation?

*How to evaluate:* Generate at 4 DRB inflow sites (Cannonsville, Pepacton, Neversink, lateral inflows) using shared bootstrap indices. Compare cross-site correlation matrix to historical. Wheeler et al. (2025) demonstrated multi-site SA generation for the Eastern Nile with explicit cross-correlation in the objective; MOEA-FIND's shared-index approach is simpler but should suffice for sites within the same basin.

*Key references:* Wheeler et al. (2025) for multi-site generation; Nowak et al. (2010) for multi-site disaggregation with K-NN conditioning; Kirsch et al. (2013) for shared-index spatial correlation.

**RQ8: Policy Vulnerability Attribution**
When combined with BART-based scenario discovery, does a MOEA-FIND ensemble reveal different or more informative vulnerability structure than conventional ensembles?

*How to evaluate:* Re-evaluate NYCOptimization Pareto policies on MOEA-FIND ensemble vs. library-subsampling baseline. Apply BART (Chipman et al. 2010) for scenario discovery. Compare identified vulnerability factors and their credible intervals. Robustness metrics should be stakeholder-specific following Hadjimichael et al. (2020) satisficing thresholds and Sunkara et al. (2023) stakeholder-contingent definitions. Structured coverage ensures BART has diverse, unbiased training data across drought characteristic space.

*Key references:* EMODPS framework for policy-scenario pairing; Chipman et al. (2010) BART; Hadjimichael et al. (2020) and Sunkara et al. (2023) for stakeholder-contingent robustness and vulnerability analysis.

---

## Figure Sequence

Figures are ordered to build the reader's understanding progressively: first the mathematical mechanism, then validation on test problems, then the hydrology application, and finally the downstream value for decision-making. Each figure directly supports one or more research questions.

**Figure 1: The Manhattan Norm Trick (Conceptual)** GENERATED (2026-04-13)
2D and 3D schematics showing: (a) the (k+1)-objective formulation with anti-ideal point, (b) the resulting hyperplane/simplex in objective space, (c) epsilon-dominance tiling on the simplex, (d) projection back to k-dimensional drought space yielding near-uniform coverage. Purpose: make the core mechanism accessible without equations. Supports RQ1.

**Figure 2: Analytic Proof-of-Concept** GENERATED (2026-04-13)
2D and 3D analytic test results. Panels: (a) 2D Pareto solutions (near-indistinguishable from LHS), (b) 3D simplex with pairwise projections showing full range coverage, (c) L2-star discrepancy and nearest-neighbor CV comparison (Pareto vs LHS vs Sobol), (d) 3D scatter showing Pareto front on hyperplane. Supports RQ1, RQ2.

**Figure 3: Dimensionality Scaling** GENERATED (2026-04-13); superseded by DD-11 dimension sweep
Coverage quality metrics (interior mass fraction, grid occupancy, orthant coverage) from K=2 to K=6 on the constrained K-ball problem. All values preliminary (EpsNSGAII stand-in; production values pending Borg HPC runs). NN_CV=0.42 and L2* values from this earlier run are superseded by the DD-11 results; do not cite them in the manuscript. See DD-10 and DD-11. Supports RQ2.

**Figure 4: Epsilon and Anti-Ideal Sensitivity**
Two-panel sensitivity analysis. (a) Coverage quality vs. epsilon value (ensemble size, spacing regularity, discrepancy). (b) Coverage quality vs. anti-ideal placement (extent of sampled drought space, fraction of epsilon-boxes filled). Demonstrates parameter selection guidance. Supports RQ6.

**Figure 5: Generator Comparison (Kirsch Bootstrap vs. Parametric)**
SSI-based drought objectives (duration, mean average severity). Panels: (a) Pareto fronts from Kirsch bootstrap, KDE, and parametric (kappa + D-vine) generators overlaid, (b) achievable drought space extent for each, (c) example synthetic hydrographs from each generator showing differences in extremes. Supports RQ3 (plausibility) and demonstrates the parametric formulation's contribution.

**Figure 6: Trace Plausibility Assessment**
Multi-panel validation of generated traces. (a) Lag-1 autocorrelation (synthetic vs. historical, by calendar month), (b) flow duration curves (historical envelope with synthetic traces overlaid), (c) seasonal cycle (monthly mean flows), (d) Hurst coefficient. Compare bootstrap (with Cholesky) vs. parametric generators. Supports RQ3.

**Figure 7: Coverage Comparison with Alternatives**
MOEA-FIND ensemble vs. (a) LHS subsampling from a large Kirsch library, (b) repeated single-objective FIND runs on a target grid, (c) Sobol nearest-neighbor selection. Metrics: L2-star discrepancy, nearest-neighbor CV, fraction of feasible drought space covered, computational cost (NFE or wall-clock). Supports RQ4.

**Figure 8: DRB Multi-Site Application**
3D drought space (frequency, duration, intensity) with MOEA-FIND Pareto front shown as a triangular simplex surface. Panels: (a) 3D scatter of ensemble members, (b) pairwise 2D projections with marginal histograms, (c) cross-site correlation matrix comparison (historical vs. synthetic). Supports RQ7.

**Figure 9: Policy Vulnerability Analysis**
Re-evaluation of NYCOptimization Pareto policies on the MOEA-FIND ensemble. (a) Robustness heatmap across policies and stakeholder-specific metrics, (b) BART partial dependence plots showing which drought characteristics drive policy failure, (c) comparison of vulnerability regions identified from MOEA-FIND ensemble vs. library-subsampling baseline. Supports RQ8.

---

## Supplemental Information (Planned)

**SI-1: Mathematical Proof of the Hyperplane Property**
Formal proof that the Pareto front of the (k+1)-objective Manhattan norm formulation lies on the (k-1)-simplex sum(J_i) = C. Cover the general case for k dimensions with conditions on anti-ideal placement.

**SI-2: SSI Calculation Details**
Full description of the SynHydro SSI pipeline: accumulation, distribution fitting per calendar month, normal transform. Drought event identification rules (onset, critical threshold, recovery hysteresis). Parameter sensitivity (SSI-1, SSI-3, SSI-6, SSI-12).

**SI-3: Generator Implementation Details**
Kirsch pipeline walkthrough (log-space deseasonalization, normal-score transform, Cholesky decomposition, Dec-Jan handling). Parametric generator details (kappa distribution fitting, D-vine copula structure, inverse Rosenblatt transform). Code availability statement.

**SI-4: Convergence Diagnostics**
Epsilon-progress and hypervolume convergence traces for all experiments. Operator selection dynamics (Borg's adaptive operator probabilities over the run). Runtime analysis (wall-clock scaling with dimensionality and NFE budget).

**SI-5: Extended Plausibility Diagnostics**
Full suite of statistical tests comparing synthetic vs. historical traces. QQ plots for monthly distributions by calendar month. Spectral analysis (power spectral density). Cross-correlation matrices for multi-site experiments (all site pairs).

**SI-6: Sensitivity to SSI Timescale**
Pareto fronts and coverage metrics for SSI-1, SSI-3, SSI-6, SSI-12 accumulation periods. Discussion of how timescale choice affects the drought characteristic space and the types of droughts the method can generate.

**SI-7: Full Pareto Front Tables**
Tabulated drought characteristics for all Pareto-optimal solutions in the DRB case study. Provided as CSV for reproducibility.
