# Critical Analysis: Things That Could Go Wrong or Need Rethinking

*This document is for honest assessment of assumptions, risks, and alternatives.*

---

## 1. The Short Trace Insight

Generating short traces (5-15 years) with 1-2 targeted drought events rather than 78-year records fundamentally changes the problem:

**What it enables:**
- Dramatically fewer DVs (60-180 monthly, or 20-60 seasonal blocks, vs. 936 for 78 years)
- Borg converges much faster in lower dimensions
- Each trace IS a drought event, not an average over many events
- Objectives become event-level (this drought's duration, this drought's intensity) rather than statistical summaries
- More interpretable: "generate a drought that lasts 8 months with peak deficit of 500 MGD"

**What it complicates:**
- Pywr-DRB needs a full simulation period. Short drought traces must be embedded in context (normal-flow years before/after). How? Options: (a) splice into historical baseline, (b) generate a "normal" baseline separately and insert the drought, (c) repeat the drought trace with normal padding.
- Initial conditions matter. Reservoir storage at drought onset affects outcomes. This is itself a dimension to explore, not to assume.
- "Drought frequency" loses meaning for a single event. The objective set must shift.

**What it suggests about the method:**
- Maybe the method is not "generate a 78-year synthetic record" but "generate a specific drought event" and the ensemble is a diverse collection of drought events, each with different characteristics
- This is closer to how stress testing actually works in practice (design specific stress scenarios)
- The embedding problem (putting the drought into a simulation context) is separate from the generation problem

---

## 2. Is the Manhattan Norm Really Necessary?

The Manhattan norm forces the Pareto front onto a hyperplane, which produces uniform coverage. But is this the only way?

**Alternative 1: Just use k objectives without the norm.**
- The Pareto front of k minimization objectives is a (k-1)-dimensional surface
- For objectives that are somewhat conflicting (e.g., high frequency vs. high duration), the Pareto front already has extent
- But without the norm, the front may not fill uniformly. It will concentrate where trade-offs are steepest.
- With the norm, uniformity is guaranteed by the simplex geometry

**Alternative 2: Use a reference-point MOEA (NSGA-III, MOEA/D) instead of the norm.**
- NSGA-III and MOEA/D use reference directions to distribute solutions evenly on the Pareto front
- Achieves similar uniformity through a different mechanism
- Does not require the extra objective
- But: these methods don't have Borg's epsilon-dominance, and the "uniformity" is along reference directions, not in the projected drought space

**Assessment:** The Manhattan norm is elegant because it converts a coverage problem into a standard Pareto optimality problem that Borg is designed to solve. The alternatives are viable but less clean. Keep the norm. The IEEE CEC (2019) literature on sampling reference points on Pareto fronts confirms that uniform distribution of solutions on irregularly shaped fronts is an active research problem (SPREAD algorithm). The Manhattan norm sidesteps this by forcing the front onto a regular hyperplane where epsilon-dominance naturally provides uniformity.

**User's position (2026-04-12):** Borg MOEA is known to be far superior to NSGA variants. Method only needs to be demonstrated with Borg. It may generalize to other MOEAs but that's not our concern for this paper.

---

## 3. Does Epsilon-Dominance Actually Produce Uniformity on Higher-Dimensional Simplices?

In 2D (line), epsilon-dominance clearly produces uniform spacing. In 3D (triangular simplex), the epsilon boxes tile the space, but the simplex intersects boxes at angles. The number of unique epsilon-boxes on a k-simplex is O(C^k / epsilon^{k-1}) where C is the simplex diameter.

**Analysis (2026-04-13, Experiment 1.2 COMPLETE):**

The 3D analytic test revealed that epsilon-dominance does tile the 2-simplex (1362 Pareto solutions satisfy sum(J_i)=C to machine precision), but coverage uniformity is moderate when measured in full DV space:

- **NN_CV (spacing regularity):** 0.42 (Pareto) vs 0.37 (LHS) vs 0.28 (Sobol) vs 0.50 (Random)
- **L2* discrepancy:** 0.038 (Pareto) vs 0.006 (LHS) vs 0.001 (Sobol) vs 0.042 (Random)

**Why the gap?** In the analytic test, all 1362 feasible points (on the [-3,3]^3 cube intersected with the hyperplane) are Pareto-optimal. There are no trade-offs; Borg's job is purely to tile the hyperplane. LHS and Sobol, by contrast, are designed to fill the full unconstrained cube uniformly. The comparison is not apples-to-apples.

**Significance for hydrology:** In real drought generation, the Pareto front will be much smaller and confined to the feasible drought region by physics, not by Borg's discretion. Epsilon-dominance will concentrate solutions where trade-offs exist and constraints bind. LHS/Sobol cannot guide this discovery a priori. The advantage of MOEA-FIND is not maximal uniformity in unconstrained space, but structured coverage of the *feasible* region that emerges from coupled generation and optimization.

**Recommendation:** Empirical validation in Phase 2 (Experiment 2.1) will test whether MOEA-FIND's coverage of the feasible Kirsch-generated drought space outperforms library-and-subsample (random Kirsch library + LHS subsampling in drought space). This is the fairest comparison.

---

## 4. Are Drought Characteristics Actually Independent?

The Manhattan norm trick works best when all objectives can vary independently. But drought characteristics are correlated:
- Longer droughts tend to be more severe (positive duration-severity correlation)
- Higher frequency means shorter inter-event times, which may limit maximum duration
- Intensity and duration may be physically linked (severe droughts exhaust baseflow, ending sooner)

**Implication:** The "achievable" region of drought characteristic space is not the full hypercube. It's a subset, possibly non-convex. The Pareto front will only span the achievable region.

**This is actually fine.** The Pareto front automatically discovers the boundary of what's achievable given the historical data and the generation mechanism. Points that are not achievable simply won't appear. The method reveals the feasible drought space as a byproduct.

**User's position (2026-04-12):** Conduct the search and see what's feasible. The correlation between drought characteristics is itself interesting to discover and report. The Pareto front's shape in drought space reveals the physically achievable trade-offs.

---

## 5. Kirsch Bootstrap vs. the Objectives: Can the Generator Actually Hit Extreme Targets?

The Kirsch bootstrap resamples blocks that exist in the historical record. The independent monthly resampling (B=1) with historical correlation imposition via Cholesky decomposition allows flexible recombination but limits extrapolation. If the worst historical drought lasted 6 months, the generator can produce droughts by recombining months from different historical droughts. But:

- Recombined sequences may have unrealistic transitions despite Cholesky enforcing correlation structure
- The historical record is finite. For DRB with 78 years of data, there are maybe 5-10 significant droughts. The combinatorial space is constrained by what months exist in the historical archive.
- Extreme targets (e.g., 3x historical duration) may be infeasible given the historical envelope

**Mitigation options:**
- Accept that the method explores the space of "rearrangements of historical drought patterns" rather than truly novel droughts
- This is actually a strength for defensibility: "all generated droughts use historical hydrological patterns as building blocks"
- Frame as: the method discovers the boundary of what's plausible given the historical record

**Two-track approach (see DD-02):** For the paper, develop both Kirsch bootstrap (respects historical envelope) and parametric CDF sampling (allows extrapolation beyond historical). Options for parametric track include: (a) fitted kappa distributions (Svensson et al. 2017) with D-vine copula for temporal dependence (Brechmann et al. 2017, Li et al. 2023), (b) KDE-smoothed empirical CDFs (mixed copula-KDE framework, 2025), (c) bootstrapped CDF uncertainty bands (BEUM approach, 2022). Phase randomization (Papalexiou 2019) is another option preserving spectral properties. Comparison between the two tracks is itself a valuable contribution.

---

## 6. What Is the Actual Computational Bottleneck?

Not Borg NFE (the generator is fast). The real bottlenecks are:

1. **Re-evaluation in Pywr-DRB** (~150s per 78-year trace). If we generate 200 drought scenarios and evaluate all Pareto policies on each, that's 200 * N_policies * 150s.
2. **Embedding short traces into simulation-ready periods.** If we generate 10-year drought events, we need to construct full simulation periods around them. This is a design problem, not a computational one.

---

## 7. Is This Actually a Better Approach Than Just Running FIND N Times?

Devil's advocate: FIND generates one targeted scenario per run. If we want 200 scenarios, run FIND 200 times with a designed grid of targets. What does MOEA-FIND add?

**MOEA-FIND advantages:**
- Discovers the feasible region automatically (FIND would fail silently at infeasible targets)
- Produces structured coverage without manual target specification
- Single run vs. 200 independent runs (though MOEA-FIND's single run is longer)
- The Pareto front has a mathematical relationship to the drought space (simplex structure)

**FIND advantages:**
- Simpler to implement
- Each scenario is guaranteed to hit its target (or fail explicitly)
- No dependency on Borg licensing
- Established in the literature

**Assessment:** MOEA-FIND is methodologically novel and produces a stronger product (structured ensemble with coverage guarantees). But for practical applications where you just need diverse droughts, FIND-with-a-grid may be "good enough." The paper needs to demonstrate when and why the MOEA-FIND approach adds value beyond FIND.

---

## 8. Publication Positioning

**Target venue (decided 2026-04-13): *Water Resources Research* (WRR).**

Rationale:
- Borgomeo et al. (2015), the most direct methodological predecessor, is a WRR paper. Reviewers familiar with search-based synthetic generation are routed there.
- The DRB multi-site case study gives the paper the applied-hydrology grounding WRR expects alongside the methods contribution.
- WRR's audience is the primary consumer of MORDM-style robustness analyses, which is where MOEA-FIND's output is intended to be used.
- Environmental Modelling & Software (where FIND, Zaniolo et al. 2023, was published) remains a credible fallback if WRR reviewers judge the scope too narrowly methodological, but the plan is written assuming WRR.

Drafting implications:
- Structure follows WRR's Methods–Results–Discussion convention; no hard word cap but target 8–12k words with 8–10 main-text display items.
- Frame the contribution as a hydrology methods advance with a case study, not as a pure optimization paper.
- Lean into the DRB case study in §7 and ensure the multi-site plausibility diagnostics are thorough (WRR reviewers will scrutinize them).

**Minimum viable paper:**
1. Method description (Manhattan norm trick, Borg coupling, bootstrap formulation)
2. Analytic proof-of-concept (2D and 3D)
3. Synthetic streamflow application (single-site, compare to FIND and library subsampling)
4. DRB case study (multi-site, show practical value)

**Stretch goals:**
- Theoretical analysis of coverage guarantees
- Scaling analysis to 5+ dimensions
- Integration with BART scenario discovery
