# notes/ INDEX — MOEA-FIND Project

*Created 2026-04-15. Living project-management document. Pointers only — edit source files, not this index.*

---

## Most Valuable / Load-Bearing Content

These items must not get lost. Each is flagged with its authoritative source.

**1. The correct L1 formulation (DD-11, design_decisions.md §DD-11)**
The method uses `f_j = D_j` (raw characteristic) and `f_{K+1} = ||D − D*||_1`. NOT `f_j = |D_j − D*_j|` — that degenerate form collapses the Pareto front to a single point. The manuscript has been corrected; any future draft work must use the DD-11 formulation.

**2. Interior coverage is empirically verified (shell_vs_interior_diagnostic.md)**
K-ball diagnostic at K = 2 through 6, 30 000 NFE each. MOEA-FIND matches uniform-in-ball on mean L1 distance at every K; interior mass fraction ≥ uniform at every K; all 2^K orthants occupied at every K. The sole degradation is grid cell coverage at K = 6 (68.4% vs 81.7% uniform), attributable to epsilon-box sub-grid clustering, not shell failure. This is the empirical foundation of the main claim.

**3. Narrative reframe required — "feasible hazard region," not "uniform coverage" (DD-10, design_decisions.md §DD-10)**
Raw coverage metrics (NN_CV = 0.42 for MOEA-FIND vs 0.37 LHS vs 0.28 Sobol) favor Sobol. This comparison is not the paper's argument. The argument is structured coverage of the *physically feasible drought hazard region*, which LHS/Sobol cannot target without a pre-generated library. The manuscript must commit to this framing; DD-10 is the authoritative statement.

**4. Terminology rules — the style_guide.md §5 banned list governs all manuscript prose**
- ALLOWED: "drought hazard space," "feasible drought hazard region," gradient boosted tree classifier
- BANNED: "outcome space," "feasibility discovery," "admissible drought characteristics," "Reed group," "BART scenario discovery" (in main text), "image of a map" (without gloss)
- These rules override planning docs and notes. style_guide.md wins on all conflicts.

**5. DD-12 is a hard gate on the empirical thesis (design_decisions.md §DD-12)**
The K-ball diagnostic does not substitute for verifying interior coverage on the Cannonsville Pareto archive (360–936 continuous DVs, non-convex feasible region). The DD-12 audit protocol (five metrics) must run after Phase beta HPC. Without it, Section 6.3 cannot be written.

**6. Apples-to-apples comparison baseline is DD-09 Option D (design_decisions.md §DD-09)**
The only valid coverage comparison is MOEA-FIND vs library-subsample in drought *characteristic* space (not in parameter space). `src/library.py` and `experiments/kirsch_ensemble/run_library_baseline.py` are ready; the 10K+ Kirsch library generation is the bottleneck. Do not compare against raw LHS/Sobol in the full bounding box.

**7. Multi-site DRB extension is explicitly deferred to follow-up (method_proposal.md sync note, publication_plan.md sync note)**
The present paper demonstrates on single-site Cannonsville only. RQ7 and Phases 4–5 are follow-up scope. Any manuscript work that reintroduces DRB multi-site as a main-text result is out of scope.

**8. Scenario discovery in §3.3 uses gradient boosted trees, not BART or PRIM (research_questions.md sync note)**
The sync note in research_questions.md and the style guide both correct this. Any manuscript text referencing BART or PRIM for the Section 3.3 demonstration is stale.

---

## File-by-File Summary

### design_decisions.md
**Status: CURRENT — primary reference for method and experiment decisions**

52 KB. Living log of 12 numbered design decisions (DD-01 through DD-12) with options, rationale, and resolution status. Authoritative on: trace length (DD-01), DV formulation (DD-02), bootstrap granularity / Kirsch B=1 (DD-03), objective set (DD-04), constraint design (DD-05), SynHydro integration via KirschBorgWrapper (DD-06), Borg interface: platypus for POC / Borg for publication (DD-07), zero-drought Pareto corner (DD-08), apples-to-apples coverage baseline in drought characteristic space (DD-09), epsilon-dominance uniformity and feasible-region reframe (DD-10), correct L1 formulation and verified interior coverage (DD-11), high-dimensional Cannonsville empirical audit protocol (DD-12).

**Overlap / source-of-truth notes:** DD-11 and shell_vs_interior_diagnostic.md contain the same coverage table — design_decisions.md is the formal record, shell_vs_interior_diagnostic.md adds context and motivation. DD-09 overlaps with the library.py / run_library_baseline.py implementation; design_decisions.md is authoritative on the decision, src/ is authoritative on the implementation.

---

### shell_vs_interior_diagnostic.md
**Status: CURRENT — empirical companion to DD-11**

11 KB. Documents the motivation, design, results table, and verdict for the K = 2–6 K-ball diagnostic. Explains the degenerate vs correct formulation, the constrained K-ball test setup, the four metrics (mean L1, interior mass fraction, orthant occupancy, grid cell coverage), the full results table, and the verdict. Also lists all output files (scripts, outputs/, figures/).

**Overlap:** The coverage table is duplicated in DD-11. This file adds the why (hostile critique history, escalation) and the diagnostic files inventory; design_decisions.md has the formal decision record. Both are current.

---

### critical_analysis.md
**Status: PARTIALLY SUPERSEDED — read with the 2026-04-14 sync note**

15 KB. Honest risk assessment covering eight topics: short vs long traces, necessity of the Manhattan norm, epsilon-dominance interior coverage, drought characteristic independence, Kirsch historical-envelope limitations, computational bottleneck, MOEA-FIND vs FIND-with-a-grid, and publication positioning.

**Superseded items:**
- §2 (Is the Manhattan Norm Necessary?) — answered by DD-11; the construction is not only necessary but verified.
- §3 (Does Epsilon-Dominance Produce Interior Coverage?) — resolved by DD-11 and shell_vs_interior_diagnostic.md.
- §5 (Kirsch Bootstrap vs. Objectives) — partially resolved; analytic case validated; hydrology case is DD-12.
- §8 (Publication Positioning) — read against style_guide.md §3.1 which defines the final section outline and venue.

**Still current:** §1 (short trace insight and embedding problem), §4 (drought characteristic independence — the Pareto front reveals the feasible region automatically), §6 (computational bottleneck: Pywr-DRB re-evaluation at ~150s/trace is the real bottleneck), §7 (MOEA-FIND vs FIND-with-a-grid — the "discover vs specify" argument remains the core contribution claim).

---

### method_proposal.md
**Status: HISTORICAL — preserved for provenance; authoritative for origin story only**

33 KB. The original design proposal (2026-04-12). Correct in spirit on the Manhattan norm (§2) but predates DD-11 and uses the wrong formulation in §2.2. Also predates the scope narrowing to single-site Cannonsville; §3.2 and Phase 4 multi-site material is deferred. The Experimental Plan (§7) is superseded by publication_plan.md Phases A–E.

**Sync note in file** identifies all superseded items; do not transcribe §3.2, Phase 3, or Phase 4 into the manuscript. The origin-story narrative for MOEA-FIND vs prior work (§1) and the architecture diagram (§3.1) remain useful as first-draft material.

**Overlap:** §5 Research Questions superseded by research_questions.md. §8 References is the most complete reference list in notes/ for the parametric track (Svensson, Brechmann, Li, Papalexiou); use as a lookup when drafting the Background section.

---

### research_questions.md
**Status: CURRENT — use alongside manuscript_main_draft.md**

13 KB. Eight research questions (RQ1–RQ8), figure sequence (Figs 1–9), and SI plan (SI-1–SI-7). The 2026-04-14 sync note in the file corrects the pre-correction state: "outcome space" → "drought hazard space"; multi-site deferred (RQ7 is a research direction, not a main-text result); §3.3 uses gradient boosted trees (not BART or PRIM); shell-vs-interior resolved as DD-11.

**Load-bearing items:** The figure sequence (Figs 1–9) is the authoritative list of planned display items. The SI structure (SI-1 through SI-7, now shifted to SI-1 through SI-8 after SI-2 was inserted for DD-11) is the reference for the supporting information draft.

---

### publication_plan.md
**Status: CURRENT with sync caveat — use Phases A–E; ignore DRB scope**

17 KB. Phased work breakdown (Phases A–E) from repo scaffolding through submission, with exit criteria for each phase and three checkpoint gates. The 2026-04-14 sync note identifies the scope narrowing (single-site Cannonsville, not DRB main-text), the five-section flat outline, the hazard-space terminology, and the gradient boosted tree correction.

**Still fully authoritative:** Phase B task list (B1–B6) and exit criteria define what is needed before HPC. Phase C task list (C0–C5) defines HPC transition steps including MM Borg setup and checkpointing. Gates 1–3 are the checkpoint criteria. Section 1.3 Known Gaps is still accurate.

**Overlap:** Phase A tasks (A1–A5) are largely complete per manuscript_main_draft.md. Section 4 Recommended Next Steps is partially stale (A1, A2 are done); the bottleneck is now the 10K-trace Kirsch library (B4) and HPC access (Phase C0).

---

### literature_review.md
**Status: CURRENT — corrected 2026-04-14; use for Background drafting**

38 KB. Structured literature review across five domains: (1) search-based generation (Borgomeo, Zaniolo/FIND, Wheeler), (2) Kirsch-Nowak block bootstrap, (3) parametric generation with copulas, (4) space-filling design and scenario subsampling (Bonham 2024 as the closest published precedent), and (5) policy re-evaluation and scenario discovery consumers.

**Key corrected items** (per 2026-04-14 sync note): Zaniolo et al. (2024) not 2023; Wheeler et al. (2025) journal is Journal of Hydrologic Engineering not JWRPM; all three predecessors (Borgomeo, FIND, Wheeler) use target-matching in parameter/DV space, not post-hoc subsampling in drought characteristic space; fabricated claims about Quinn, Hadjimichael, Herman using hazard-space subsampling were removed.

**Overlap with method_proposal.md §8:** method_proposal.md §8 contains the most complete parametric-track reference list. For Background §2 (predecessor characterization), literature_review.md §1 is the source of truth.

---

### design_decisions.md / notes/literature/
**Status: CURRENT — four micro-summaries of core papers**

4 files, ~6 KB total: borgomeo_2015.md, parametric_generation.md, wheeler_2025.md, zaniolo_2024_FIND.md. Each is a short (1–2 page) focused summary of one key paper. Useful for quick lookups when drafting.

---

### style_guide.md
**Status: CURRENT — MANDATORY for all manuscript prose; highest priority**

29 KB. Nine explicit user rules (§2.1–2.8) governing bold/italic, code formatting, informal language, technical term definition, prose vs lists, section headers, narrative sequencing, and overstatement. Literature-grounded rules (§3.1–3.10) covering section headers, paragraph structure, math notation, lists, em-dashes, abstract structure, captions, and transitions. Domain vocabulary reference (§4) with allowed/disallowed terms. Anti-pattern checklist (§5.1–5.17+) with the banned-terms list for manuscript prose.

**This file wins on all rule conflicts.** Any planning note, design decision, or prior draft that contains banned terminology (§5.11–5.17: "outcome space," "feasibility discovery," "admissible drought characteristics," "Reed group," "BART scenario discovery") must not have that language transcribed into the manuscript. The guide specifies replacement phrases.

---

## Cross-Reference Map

| Topic | Primary source | Secondary / corroborating |
|---|---|---|
| Correct L1 formulation | design_decisions.md §DD-11 | shell_vs_interior_diagnostic.md §"underlying issue" |
| Interior coverage evidence | shell_vs_interior_diagnostic.md §Results | design_decisions.md §DD-11, critical_analysis.md §3 |
| Coverage uniformity vs feasible-region framing | design_decisions.md §DD-10 | publication_plan.md §1.3 gap 2 |
| Apples-to-apples comparison baseline | design_decisions.md §DD-09 | publication_plan.md §B4 |
| Trace length / event vs aggregate | design_decisions.md §DD-01 | critical_analysis.md §1, research_questions.md RQ5 |
| DV formulation (index vs residual injection) | design_decisions.md §DD-06 | method_proposal.md §3.0 |
| Constraint design | design_decisions.md §DD-05 | method_proposal.md §6.3 |
| Borg vs platypus | design_decisions.md §DD-07 | publication_plan.md §B2 |
| Predecessors (Borgomeo, FIND, Wheeler) | literature_review.md §1 | method_proposal.md §8 |
| Multi-site deferred | method_proposal.md sync note | publication_plan.md sync note, research_questions.md RQ7 |
| Scenario discovery = gradient boosted trees | research_questions.md sync note | style_guide.md §5 |
| Manuscript prose rules | style_guide.md §2–5 | (overrides everything) |
| Phased work plan / exit criteria | publication_plan.md §2 | — |
| DD-12 empirical audit protocol | design_decisions.md §DD-12 | — |

---

## Flagged Contradictions and Stale Content

These are flags only — source files are not edited.

**FLAG 1 — Algorithm identity inconsistency (TREV-DECISION)**
`design_decisions.md §DD-07` defers the Borg vs platypus question ("use platypus for POC, Borg for publication"). The analytic experiments currently use EpsNSGAII (platypus). The manuscript names Borg throughout. `methods_audit.md` Discrepancy 5 documents this as submission-blocking. The section outline in `style_guide.md §3.1` names a nine-section structure that still includes multi-site DRB (sections 7 and 8); this section list has not been updated to match the scope narrowing in publication_plan.md. Both require Trev decision before manuscript lock.

**FLAG 2 — Archive size mismatch across documents**
`CLAUDE.md` (project root) states 1362 Pareto solutions; shell_vs_interior_diagnostic.md and DD-11 report 6158 at K = 3 (30 000 NFE) and 682 at K = 3 (2 000 NFE). method_proposal.md §7 Experiment 1.2 cites 1362 from a prior run. These are from different NFE budgets, not errors, but any manuscript statement must specify the run configuration. `methods_audit.md` Discrepancy 7 flags this.

**FLAG 3 — Constraint tolerance mismatch (TREV-DECISION)**
`methods_audit.md` Discrepancies 1–4 document manuscript-vs-code mismatches: lag-1 AC tolerance (manuscript 0.05 vs code 0.30), non-drought mean (manuscript 15%/non-drought vs code 50%/all-flows), seasonal cycle constraint (in code, absent from manuscript), annual mean/CV tolerance (50% in code, undocumented in manuscript). These are TREV-DECISION items that must be resolved before the manuscript Methods section can be finalized. design_decisions.md §DD-05 has the design intent but not the implemented values.

**FLAG 4 — Section outline version drift**
`style_guide.md §3.1` lists a nine-section outline (Introduction, Background, Methods, Study Area, Analytic Validation, Single-Site, Multi-Site, Discussion, Conclusions) that still includes multi-site DRB (sections 7–8). `publication_plan.md` and all sync notes state multi-site is deferred. The authoritative outline is in `manuscript_main_draft.md`, which should be treated as the master. style_guide.md §3.1 is descriptive, not prescriptive, on section count.

**FLAG 5 — "Hazard-space is emergent from parameter-space" framing absent from notes**
`research_questions.md`, `method_proposal.md`, and `critical_analysis.md` all frame the novelty as sampling in "drought outcome/characteristic space." Neither explicitly states the key positioning nuance that hazard-space coverage *emerges* from the search over parameter space — it is not directly sampled in parameter space (contrast with Hadjimichael 2020, Bonham 2024). The framing in `manuscript_main_draft.md` Introduction carries this nuance; the notes do not. When revising any notes section for manuscript prose, explicitly add this distinction.

**FLAG 6 — Moallemi Fork 3.2.2 / Herman Axis II vocabulary not in notes**
The Herman Axis II / Moallemi Fork 3.2.2 framing (referenced in the CLAUDE.md framing for the framing_anchor_synthesis.md) does not appear in any notes file. If this framing is load-bearing for the Discussion, it lives only in framing_anchor_synthesis.md (a root-level manuscript file) and has not been synthesized into the notes/. Consider adding a pointer to framing_anchor_synthesis.md in this INDEX when the Discussion is drafted.

---

## Retired Files

*All files evaluated 2026-04-15. No files retired.*

Every file in notes/ was evaluated against the criteria: (a) classified as "historical" or "superseded" in the summary above, and (b) all load-bearing content preserved in INDEX.md or a current source-of-truth file.

**method_proposal.md** — Classified HISTORICAL but retained. Unique content not preserved elsewhere: §3.1 architecture diagram (ASCII art of the Borg ↔ KirschBorgWrapper ↔ SynHydro loop), §8 parametric-track reference list (Svensson 2017, Brechmann 2017, Li 2023, Papalexiou 2019 and several others not fully reproduced in literature_review.md). Deleting would lose these references. Retained per "when in doubt, keep."

**critical_analysis.md** — Classified PARTIALLY SUPERSEDED but retained. §1 (short-trace insight and embedding problem), §4 (drought characteristic independence analysis), §6 (computational bottleneck), §7 (MOEA-FIND vs FIND-with-a-grid argument) are current and unique. Deleting would lose the §7 argument which is the core claim for the Discussion. Retained.
