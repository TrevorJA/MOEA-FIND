# notes/ INDEX — MOEA-FIND Project

*Created 2026-04-15. Revised 2026-04-15 after Phase C clarification propagation.*
*Living project-management document. Pointers only — edit source files, not this index.*

---

## Most Valuable / Load-Bearing Content

**1. The correct L1 formulation (DD-11, design_decisions.md §DD-11)**
The method uses `f_j = D_j` (raw characteristic) and `f_{K+1} = ||D − D*||_1`. NOT `f_j = |D_j − D*_j|` — that degenerate form collapses the Pareto front to a single point. The manuscript has been corrected; any future draft work must use the DD-11 formulation.

**2. Interior coverage is empirically verified — preliminary (shell_vs_interior_diagnostic.md)**
K-ball diagnostic at K = 2 through 6, 30 000 NFE each, using EpsNSGAII as a local stand-in. Results are preliminary pending Borg MOEA production runs on HPC. The coverage findings are valid for the construction because both algorithms use the same epsilon-dominance archive mechanism. MOEA-FIND matches uniform-in-ball on mean L1 distance at every K; interior mass fraction ≥ uniform at every K; all 2^K orthants occupied at every K. Grid cell coverage at K = 6 is 68.4% vs 81.7% uniform — attributable to epsilon-box sub-grid clustering, not shell failure.

**3. Narrative framing — "feasible hazard region coverage," not "uniform coverage" (DD-10)**
Raw coverage metrics (older exploratory run, superseded by DD-11) favored Sobol. The paper's argument is structured coverage of the *physically feasible drought hazard region*, which LHS/Sobol cannot target without a pre-generated library. The manuscript commits to this framing; DD-10 is the authoritative statement. All preliminary numbers from DD-11 and the analytic diagnostic will be confirmed with Borg production runs before submission.

**4. Terminology rules — style_guide.md §5 banned list governs all manuscript prose**
- ALLOWED: "drought hazard space," "feasible drought hazard region," gradient boosted tree classifier
- BANNED: "outcome space," "feasibility discovery," "admissible drought characteristics," "Reed group," "BART scenario discovery" (main text), "image of a map" (without gloss)
- style_guide.md wins on all conflicts.

**5. DD-12 is a hard gate on the empirical thesis (design_decisions.md §DD-12)**
The K-ball diagnostic does not substitute for verifying interior coverage on the Cannonsville Pareto archive (360–936 continuous DVs, non-convex feasible region). DD-12 audit must run after HPC Phase beta. Without it, Section 3.2 cannot be written.

**6. Apples-to-apples comparison baseline is DD-09 Option D (design_decisions.md §DD-09)**
The only valid coverage comparison is MOEA-FIND vs library-subsample in drought *characteristic* space (not parameter space). `src/library.py` is ready; the 10K+ Kirsch library is the bottleneck. Do not compare against raw LHS/Sobol in the full bounding box.

**7. Multi-site DRB extension is explicitly deferred to follow-up**
The present paper demonstrates on single-site Cannonsville only. Multi-site is native to the Kirsch shared-indices mechanism and worth noting in Methods, but is not a main-text result. RQ7 and Phases 4–5 are follow-up scope.

**8. Herman Axis II + Moallemi Fork 3.2.2 vocabulary is load-bearing for framing**
MOEA-FIND operates at Axis II (States of the World) in Herman et al. (2015) and Stage II fork 3.2.2 (Framing future scenarios) in Moallemi et al. (2020b). This vocabulary has been propagated to: manuscript_main_draft.md §4, introduction_revised.md MOEA-FIND positioning paragraph, reviewer_critiques.md Critique 9, critical_analysis.md sync update. Authoritative reference: `framing_anchor_synthesis.md`.

**9. Constraints verify statistical plausibility — not drought identification**
The two plausibility constraints verify that each synthetic trace lies within the envelope of traces the Kirsch generator could produce naturally given an infinite ensemble. They do not identify drought events. The exact tolerance values are code alignment backlog items (see `methods_audit.md` Items 1 and 2). The manuscript is the specification; the code aligns to it.

**10. Borg is the production algorithm; EpsNSGAII is the local stand-in**
All analytic benchmark experiments were run with EpsNSGAII (platypus) for rapid local development. Borg MOEA is the production algorithm. The theoretical argument applies to any epsilon-dominance archive, so analytic results are valid for the construction. All numerical values from analytic runs are preliminary pending Borg production runs on HPC.

---

## File-by-File Summary

### design_decisions.md
**Status: CURRENT — primary reference for method and experiment decisions**

52 KB. Living log of 12 numbered design decisions (DD-01 through DD-12) with options, rationale, and resolution status. Authoritative on: trace length (DD-01), DV formulation (DD-02), bootstrap granularity / Kirsch B=1 (DD-03), objective set (DD-04), constraint design as plausibility envelope verification (DD-05, resolved per `manuscript_main_draft.md §2.2.4`), SynHydro integration via KirschBorgWrapper (DD-06), Borg interface: EpsNSGAII as local stand-in / Borg as production algorithm (DD-07, updated 2026-04-15), zero-drought Pareto corner (DD-08), apples-to-apples coverage baseline in drought characteristic space (DD-09), epsilon-dominance uniformity and feasible-region reframe (DD-10), correct L1 formulation and verified interior coverage — preliminary (DD-11), high-dimensional Cannonsville empirical audit protocol (DD-12).

---

### shell_vs_interior_diagnostic.md
**Status: CURRENT — empirical companion to DD-11; all values preliminary**

11 KB. Documents the motivation, design, results table, and verdict for the K = 2–6 K-ball diagnostic using EpsNSGAII. Results are preliminary pending Borg production runs. The verdict (interior coverage empirically confirmed, shell-only and orthant-collapse failure modes refuted) is valid for the construction because both algorithms use the same epsilon-dominance archive. Output files listed.

---

### critical_analysis.md
**Status: PARTIALLY SUPERSEDED — read with the 2026-04-14 sync update**

15 KB. Honest risk assessment covering eight topics. Superseded: §2 (Manhattan norm necessity — resolved by DD-11), §3 (epsilon-dominance interior coverage — resolved by DD-11/diagnostic, values preliminary), §5 (Kirsch bootstrap — analytic case validated; hydrology case is DD-12). Still current: §1 (short-trace embedding), §4 (drought characteristic independence), §6 (computational bottleneck), §7 (MOEA-FIND vs FIND — core contribution claim).

Sync update now includes Herman Axis II + Moallemi Fork 3.2.2 framing references.

---

### method_proposal.md
**Status: ARCHIVED 2026-04-15 — preserved for provenance only**

33 KB. Original design proposal (2026-04-12). Predates DD-11 and uses the wrong formulation in §2.2. KEEP-HISTORICAL banner in file. Unique retained content: §3.1 architecture diagram (ASCII art of the Borg ↔ KirschBorgWrapper ↔ SynHydro loop), §8 parametric-track reference list (Svensson 2017, Brechmann 2017, Li 2023, Papalexiou 2019). Do not transcribe §3.2, Phase 3, Phase 4 multi-site content into the manuscript.

---

### research_questions.md
**Status: CURRENT — use alongside manuscript_main_draft.md**

13 KB. Eight research questions (RQ1–RQ8), figure sequence (Figs 1–9), SI plan (SI-1–SI-7). Sync note corrects pre-consolidation state. Figure 3 description updated 2026-04-15 to note EpsNSGAII stand-in and superseded numbers. Load-bearing: the planned figure sequence and SI structure.

---

### publication_plan.md
**Status: CURRENT — use Phases A–E; ignore DRB scope**

17 KB. Phased work breakdown (Phases A–E) with exit criteria and three checkpoint gates. §1.3 Known Gaps item 5 updated 2026-04-15: EpsNSGAII is confirmed local stand-in; Borg installation on HPC is the production prerequisite. Phases B and C task lists are fully authoritative. Section 4 Recommended Next Steps partially stale (A1, A2 done; bottleneck is now 10K-trace Kirsch library and HPC access).

---

### literature_review.md
**Status: CURRENT — corrected 2026-04-14; use for Background drafting**

38 KB. Structured review across five domains: search-based generation (Borgomeo, Zaniolo/FIND, Wheeler), Kirsch-Nowak bootstrap, parametric generation, space-filling design and scenario subsampling (Bonham 2024), policy re-evaluation and scenario discovery. All three predecessor methods (Borgomeo, FIND, Wheeler) confirmed as target-matching in parameter/DV space; fabricated hazard-space subsampling claims removed.

---

### notes/literature/ (four micro-summaries)
**Status: CURRENT — for quick lookups when drafting**

borgomeo_2015.md, parametric_generation.md, wheeler_2025.md, zaniolo_2024_FIND.md. Short focused summaries for Background drafting.

---

### style_guide.md
**Status: CURRENT — MANDATORY for all manuscript prose; highest priority**

29 KB. Nine user rules (§2.1–2.8) and literature-grounded rules (§3.1–3.10). Domain vocabulary reference (§4). Anti-pattern checklist (§5) with banned terms. This file wins on all rule conflicts.

---

## Cross-Reference Map

| Topic | Primary source | Secondary |
|---|---|---|
| Correct L1 formulation | design_decisions.md §DD-11 | shell_vs_interior_diagnostic.md §"underlying issue" |
| Interior coverage evidence (preliminary) | shell_vs_interior_diagnostic.md §Results | design_decisions.md §DD-11, critical_analysis.md §3 |
| Feasible-region framing | design_decisions.md §DD-10 | publication_plan.md §1.3 gap 2 |
| Apples-to-apples comparison | design_decisions.md §DD-09 | publication_plan.md §B4 |
| Trace length / event vs aggregate | design_decisions.md §DD-01 | critical_analysis.md §1, research_questions.md RQ5 |
| DV formulation (index vs residual) | design_decisions.md §DD-06 | method_proposal.md §3.0 |
| Constraint plausibility design | design_decisions.md §DD-05 | methods_audit.md Items 1–4 |
| Code alignment backlog | methods_audit.md | design_decisions.md §DD-05, §DD-07 |
| Borg as production algorithm | design_decisions.md §DD-07 | publication_plan.md §1.3 item 5 |
| Herman Axis II + Moallemi Fork 3.2.2 | framing_anchor_synthesis.md | manuscript_main_draft.md §4, introduction_revised.md |
| Predecessors (Borgomeo, FIND, Wheeler) | literature_review.md §1 | method_proposal.md §8 |
| Multi-site deferred | publication_plan.md sync note | research_questions.md RQ7 |
| Scenario discovery = gradient boosted trees | research_questions.md sync note | style_guide.md §5 |
| Manuscript prose rules | style_guide.md §2–5 | (overrides everything) |
| Phased work plan / exit criteria | publication_plan.md §2 | — |
| DD-12 empirical audit protocol | design_decisions.md §DD-12 | — |

---

## Active Flags

**FLAG A — Archive size mismatch across documents**
CLAUDE.md (project root) states 1362 Pareto solutions; DD-11 and shell_vs_interior_diagnostic.md report 6158 at K = 3 (30 000 NFE) and 682 at K = 3 (2 000 NFE) — all preliminary (EpsNSGAII stand-in). The 1362 figure derives from an early exploratory run. All manuscript statements must cite run configuration explicitly; the DD-11 table is the reference. `methods_audit.md` Item 7 flags the CLAUDE.md update.

**FLAG B — Section outline version drift**
style_guide.md §3.1 lists a nine-section outline that still includes multi-site DRB sections. The authoritative outline is in `manuscript_main_draft.md`. style_guide.md §3.1 is descriptive, not prescriptive, on section count.

**FLAG C — "Hazard-space is emergent from parameter-space" nuance**
Several notes frame novelty as "sampling in drought outcome space" without stating explicitly that the hazard space is emergent from the parameter space — MOEA-FIND searches parameter space to achieve hazard coverage, not the reverse. `framing_anchor_synthesis.md §4` and `manuscript_main_draft.md §4` carry this nuance. When migrating any notes text into the manuscript, add this distinction.

---

## Resolved (formerly flagged, now dissolved)

**Resolved 1 — Algorithm identity (was FLAG 1)**
EpsNSGAII is the confirmed local stand-in; Borg MOEA is the production algorithm. No discrepancy. DD-07 updated 2026-04-15. methods_audit.md reframed accordingly (Item 5 is a note, not a discrepancy).

**Resolved 2 — Constraint tolerance alignment (was FLAG 3)**
Manuscript is the specification; code is the implementation. Constraint tolerance values are code alignment backlog items (methods_audit.md Items 1–4), not TREV-DECISIONS about changing the manuscript. DD-05 updated 2026-04-15 with plausibility envelope framing.

**Resolved 3 — Herman Axis II + Moallemi Fork 3.2.2 vocabulary (was FLAG 6)**
Propagated to: manuscript_main_draft.md §4, introduction_revised.md (MOEA-FIND positioning paragraph), reviewer_critiques.md (Critique 9 defensible response), critical_analysis.md (sync update), framing_anchor_synthesis.md (source of truth, unchanged).

---

## Retired Files

**method_proposal.md** — ARCHIVED 2026-04-15. KEEP-HISTORICAL banner in file. Retained for §3.1 architecture diagram and §8 parametric reference list.

*All other files in notes/ evaluated 2026-04-15 — none retired.*
