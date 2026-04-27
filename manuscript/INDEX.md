# manuscript/ INDEX

*Revised 2026-04-27 — function-based subfolders. See `CONVENTIONS.md` for
folder rules and maintenance practices. Pointers only; edit source files,
not this index.*

---

## Master doc index

| File | What it governs |
|------|-----------------|
| `drafts/manuscript_main_draft.md` | Full prose draft — never delete. |
| `drafts/supporting_info_draft.md` | Supporting information draft — never delete. |
| `governance/style_guide.md` | Prose rules, vocabulary, banned terms. Mandatory; wins all conflicts. |
| `governance/design_decisions.md` | DD-01..DD-N: all method and experiment decisions (living doc). |
| `reference/framing_anchor.md` | Herman Axis II + Moallemi fork 3.2.2 framing; DMDU vocabulary. |
| `reference/literature_review.md` | Five-domain structured literature review. |
| `reference/reviewer_defenses.md` | Anticipated reviewer critiques and defensible responses (HC-1..HC-5). |
| `planning/publication_plan.md` | Phase tracker (A–E), exit criteria, research questions (RQ1..RQ6). |
| `planning/code_alignment_backlog.md` | Manuscript→code alignment items. |
| `planning/experiment_plan.md` | Section 3.2 / 3.3 experiment design and pipeline pointer. |
| `planning/code_state.md` | Snapshot of `src/` and `workflows/`; refresh after HPC pulls. |
| `evidence/shell_vs_interior_diagnostic.md` | K-ball K=2..6 dimension sweep; DD-11 evidence. |

---

## Cross-reference: where to find things

| Topic | Primary | Secondary |
|-------|---------|-----------|
| Correct L1 formulation | `design_decisions.md §DD-11` | `shell_vs_interior_diagnostic.md` |
| Interior coverage evidence (preliminary) | `shell_vs_interior_diagnostic.md` | `design_decisions.md §DD-11` |
| Feasible-region framing | `design_decisions.md §DD-10` | `publication_plan.md §1` |
| Apples-to-apples comparison baseline | `design_decisions.md §DD-09` | `publication_plan.md Phase B` |
| Production constraint choice (AD) | `design_decisions.md §DD-14` | `experiment_plan.md`, exp13/14 SI |
| Constraint history (superseded) | `design_decisions.md §DD-05` | (points to DD-14) |
| Code alignment backlog | `code_alignment_backlog.md` | `code_state.md` |
| Borg as production algorithm | `design_decisions.md §DD-07` | `code_state.md` |
| Herman Axis II + Moallemi fork 3.2.2 | `framing_anchor.md` | `manuscript_main_draft.md §4` |
| Predecessors (Borgomeo, FIND, Wheeler) | `literature_review.md §1` | `literature/notes/` |
| §3.3 redesign (Pywr-DRB failure labels) | `experiment_plan.md` Part A | `reviewer_defenses.md` HC-2 |
| HPC code state | `code_state.md` | `code_alignment_backlog.md` Item 10 |
| DD-12 empirical audit protocol | `design_decisions.md §DD-12` | — |
| Multi-site deferred | `publication_plan.md §4` | — |
| Manuscript prose rules | `style_guide.md §2–5` | (overrides everything) |

---

## Active flag

**Hazard-from-parameter nuance.** MOEA-FIND searches parameter space to
achieve hazard-space coverage, not the reverse. When migrating any notes
text into the manuscript, preserve this distinction. See
`framing_anchor.md §4` and `manuscript_main_draft.md §4`.

---

## Other locations

`literature/notes/` — short paper summaries: `borgomeo_2015.md`,
`parametric_generation.md`, `wheeler_2025.md`, `zaniolo_2024_FIND.md`.

`scratch/` — figure-1 critique log, figure-1 spec, intro framing moves.
Ephemeral; safe to delete once merged into a master doc.

`archive/` — KEEP-HISTORICAL files only. Currently holds
`critical_analysis_v1.md` (the original critical analysis; superseded
content has been folded into `reviewer_defenses.md` and `design_decisions.md`).

`figures/` — figure scripts and generated outputs.
