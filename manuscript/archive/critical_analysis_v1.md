# Critical Analysis (v1 — KEEP-HISTORICAL)

> **KEEP-HISTORICAL — archived 2026-04-27.** This document was the
> early honest-assessment write-up. The active concerns it raised have
> migrated as follows:
>
> - §2 (Manhattan norm necessity) and §3 (interior coverage at higher
>   K) are resolved by `governance/design_decisions.md` §DD-11 and
>   `evidence/shell_vs_interior_diagnostic.md`.
> - §1 (short-trace insight) is folded into `governance/design_decisions.md`
>   §DD-01.
> - §5 (Kirsch envelope) and §7 (better than running FIND N times) are
>   addressed by `reference/reviewer_defenses.md` Critiques 7, 12, and
>   18.
> - §6 (computational bottleneck) is addressed by
>   `reference/reviewer_defenses.md` HC-5.
> - §8 (publication positioning) is superseded by the finalised framing
>   in `reference/framing_anchor.md` and `drafts/manuscript_main_draft.md`.
>
> The original prose is retained below for provenance only. Do not edit.

*This document is for honest assessment of assumptions, risks, and alternatives.*

> **2026-04-14 sync update.** Two of the critical questions in this
> document have been answered empirically since it was written and
> should be read with the resolutions in mind:
>
> - Section 2 "Is the Manhattan Norm Really Necessary?" and Section
>   3 "Does Epsilon-Dominance Actually Produce Interior Coverage at
>   Higher K?" are resolved by DD-11 and the $K = 2$ through $K = 6$
>   shell-versus-interior dimension sweep documented in
>   `shell_vs_interior_diagnostic.md`. The construction produces
>   interior-filling coverage and populates every signed orthant on
>   the constrained $K$-ball benchmark. A residual empirical audit
>   on the high-dimensional Cannonsville case is scoped as DD-12
>   and is blocked on HPC Phase beta.
> - Section 5 "Kirsch Bootstrap vs. the Objectives" is partially
>   resolved: the construction is mathematically sound on the
>   analytic benchmark, and the realistic-hydrology question is
>   exactly what DD-12 audits when Phase beta lands.
>
> Section 8 "Publication Positioning" should be read against the
> finalised manuscript framing, which is that MOEA-FIND samples
> directly in drought hazard space rather than in parameter or
> decision-variable space (contrast with Bonham 2024 input-space
> subsampling, Hadjimichael 2020 LHS over HMM hyperparameters, and
> the broader exploratory modelling tradition reviewed by Moallemi
> 2020). In the Herman et al. (2015) four-axis taxonomy, MOEA-FIND
> operates within Axis II (States of the World), introducing directed
> search in hazard-outcome space as a new mode of scenario generation.
> In the Moallemi et al. (2020b) decision support framework, it
> operates at Stage II fork 3.2.2 (Framing future scenarios). The
> scenario discovery demonstration in Section 3.3 of the manuscript
> uses a gradient boosted tree classifier, not BART or PRIM.
> Terminology rules added to `style_guide.md` sections 5.11 through
> 5.17 apply to any prose that migrates from this file into the
> manuscript.

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

The Manhattan norm is the single construction that makes the method
work. Without it, `K` minimisation objectives with no auxiliary are
either trivially resolved (the Pareto front collapses to the
ideal point if all objectives are minimised toward benign values)
or highly irregular in shape. The auxiliary objective places every
feasible point on a known `K`-dimensional affine subset of the
`K + 1` objective space, makes every pair of feasible points
mutually non-dominated under minimisation, and makes the projection
back to drought characteristic space a bijection that preserves
dimension and orthogonality. None of these three properties is
available without the auxiliary objective.

**Status (updated 2026-04-14, DD-11).** An opus method critique
flagged an earlier version of manuscript section 3.2 as describing
a degenerate formulation in which `f_j = |D_j - D*_j|`, under which
the Pareto front collapses to the anti-ideal and the method fails.
The critique was correct about the prose but inaccurate about the
code: the implementation uses `f_j = D_j` with `f_{K+1} = ||D - D*||_1`,
which is the formulation described above and does work. The
manuscript section 3.2 and SI-1 were rewritten to reflect the
actual formulation, and an empirical diagnostic was added
(SI-2, see `shell_vs_interior_diagnostic.md`) to confirm
interior-filling coverage on constrained `K`-ball test problems
from `K = 2` through `K = 6`.

**Alternatives still viable but not adopted.** Reference-direction
methods (NSGA-III, Deb and Jain 2014) and decomposition methods
(MOEA/D, Zhang and Li 2007) distribute solutions on a Pareto front
through user-supplied directions or scalarising subproblems. They
can be made to deliver similar coverage, but they do not take
advantage of a pre-aligned feasible image and they require the
user to know the Pareto front geometry in advance. The Manhattan
norm construction pre-aligns the feasible image to a known
`K`-dimensional affine subset, so Borg's epsilon-box archive tiles
it without any front-shape knowledge. This is why the present work
restricts itself to Borg.

---

## 3. Does Epsilon-Dominance Actually Produce Interior Coverage at Higher `K`?

The original form of this question was whether Borg's epsilon-box
archive produces uniform coverage of the `(K - 1)`-simplex that the
earlier manuscript prose described. That framing was based on a
formulation the code does not execute and is no longer the right
question. The question actually addressed by the existing experiments
is whether the `K`-dimensional feasible image in the correct
formulation is tiled interior-fillingly by Borg's archive.

**Analysis (2026-04-14, DD-11 dimension sweep complete; preliminary — EpsNSGAII stand-in, pending Borg production runs).** A constrained `K`-ball diagnostic (`scripts/diag_shell_vs_interior.py`) was run at `K = 2, 3, 4, 5, 6` with 30 000 NFE per `K`, using the actual implemented formulation and comparing MOEA-FIND to uniform-in-ball, LHS-in-ball, and Sobol-in-ball reference samples at matched archive size. The full table is at `outputs/diag_shell_vs_interior/sweep_table.md` and in SI-2.4. The key metrics are (all values preliminary pending Borg production runs):

- **Mean Manhattan distance from the anti-ideal.** MOEA-FIND tracks
  the uniform-in-ball reference within sampling error at every
  `K`. No shell-only concentration at any tested `K`.
- **Orthant occupancy.** MOEA-FIND populates every one of the
  `2^K` signed orthants at every `K`, including 64 of 64 at `K = 6`.
- **Interior mass fraction.** MOEA-FIND matches the uniform-in-ball
  reference at `K in {2, 3}` and exceeds it at `K in {4, 5, 6}`.
  At `K = 6` the interior mass fraction is 0.648 for MOEA-FIND
  versus 0.540 for uniform.
- **Grid cell coverage.** Matches uniform at `K <= 5` (all at
  0.97 to 1.00). At `K = 6` with 6 bins per axis, MOEA-FIND
  covers 0.684 of the 3072 feasible cells while uniform covers
  0.817. The deviation is an epsilon-dominance sub-grid
  clustering effect in the interior of the feasible region, not
  a shell failure, and it can be mitigated with smaller epsilon
  or more NFE.

**Interpretation.** The interior-filling coverage claim is
empirically verified on the analytic `K`-ball problem at every
dimensionality from `K = 2` through `K = 6`. The fine-scale
clustering at `K = 6` is a real observation but does not affect
interior coverage or orthant occupancy. Main section 3.2 and SI-1
describe the construction; SI-2.4 reports the empirical sweep;
DD-11 summarises the diagnostic findings.

**Remaining empirical risk.** The diagnostic is analytic and has a
low-dimensional decision space (2 to 6 decision variables). The
real Kirsch single-site case has 360 to 936 continuous decision
variables and a non-convex feasible drought characteristic region.
The analogous empirical audit on the Kirsch case is part of main
section 6.3 and is a hard gate on the empirical thesis. An
orthant-occupancy diagnostic has been added to the Cannonsville
analysis pipeline to detect the same failure modes that SI-2
rules out on the analytic problem.

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
- Environmental Modelling & Software (where FIND, Zaniolo et al. 2024, was published) remains a credible fallback if WRR reviewers judge the scope too narrowly methodological, but the plan is written assuming WRR.

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
