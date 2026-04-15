# manuscript/ INDEX

*Revised 2026-04-15 — reorganized structure. See `ORGANIZATION.md` for folder rules.*
*Pointers only. Edit source files, not this index.*

---

## Non-Negotiables (read before touching anything)

1. **Correct L1 formulation (design_decisions.md §DD-11):** `f_j = D_j`, `f_{K+1} = ||D − D*||_1`. NOT `|D_j − D*_j|`.
2. **style_guide.md wins on all prose conflicts.** Banned terms: "outcome space," "feasibility discovery," "Reed group," "BART scenario discovery" (main text).
3. **Manuscript is the specification. Code aligns to it** — not the reverse. See `code_alignment_backlog.md`.
4. **All analytic results are preliminary** (EpsNSGAII local stand-in; pending Borg HPC runs). Say so whenever citing numbers.
5. **DD-12 is a hard gate** — Section 3.2 cannot be written until the Cannonsville empirical audit runs.

---

## Master Doc Index

| File | What it governs |
|---|---|
| `manuscript_main_draft.md` | Full prose draft — never delete |
| `supporting_info_draft.md` | Supplemental information draft — never delete |
| `style_guide.md` | Prose rules, vocabulary, banned terms (MANDATORY) |
| `design_decisions.md` | DD-01 through DD-12: all method/experiment decisions |
| `framing_anchor.md` | Herman Axis II + Moallemi Fork 3.2.2 framing; source of truth for DMDU vocabulary |
| `reviewer_defenses.md` | Critique list with defensible responses |
| `publication_plan.md` | Phased work plan (Phases A–E), exit criteria, known gaps |
| `code_alignment_backlog.md` | Manuscript→code alignment items; not a discrepancy log |
| `constraints_spec.md` | Constraint design: what Kirsch preserves, gap analysis, calibration recommendations |
| `critical_analysis.md` | Risk assessment; partially superseded — read with sync update notes |
| `research_questions.md` | RQ1–RQ8, figure sequence (Figs 1–9), SI plan |
| `shell_vs_interior_diagnostic.md` | K-ball diagnostic K=2–6 — all values preliminary |
| `literature_review.md` | Five-domain structured literature review |
| `experiment_plan.md` | HPC deployment tracker; §3.3 redesign (Pywr-DRB failure labels, Hashimoto metrics) |
| `ORGANIZATION.md` | Folder layout, file roles, sweep/archive policy |

---

## Cross-Reference: Where to Find Things

| Topic | Primary | Secondary |
|---|---|---|
| Correct L1 formulation | `design_decisions.md §DD-11` | `shell_vs_interior_diagnostic.md §underlying issue` |
| Interior coverage evidence (preliminary) | `shell_vs_interior_diagnostic.md §Results` | `design_decisions.md §DD-11` |
| Feasible-region framing | `design_decisions.md §DD-10` | `publication_plan.md §1.3 gap 2` |
| Apples-to-apples comparison baseline | `design_decisions.md §DD-09` | `publication_plan.md §B4` |
| Constraint plausibility design | `design_decisions.md §DD-05` | `constraints_spec.md` |
| Code alignment backlog | `code_alignment_backlog.md` | `design_decisions.md §DD-05, §DD-07` |
| Borg as production algorithm | `design_decisions.md §DD-07` | `publication_plan.md §1.3 item 5` |
| Herman Axis II + Moallemi Fork 3.2.2 | `framing_anchor.md` | `manuscript_main_draft.md §4` |
| Predecessors (Borgomeo, FIND, Wheeler) | `literature_review.md §1` | `archive/method_proposal.md §8` |
| §3.3 redesign (Pywr-DRB failure labels) | `experiment_plan.md §Part A` | `reviewer_defenses.md HC-2` |
| HPC deployment TODOs | `experiment_plan.md §Part B` | — |
| DD-12 empirical audit protocol | `design_decisions.md §DD-12` | — |
| Multi-site deferred | `publication_plan.md sync note` | `research_questions.md RQ7` |
| Manuscript prose rules | `style_guide.md §2–5` | (overrides everything) |

---

## Active Flags

**FLAG A — Archive size mismatch across documents.**
CLAUDE.md states 1362 Pareto solutions; `design_decisions.md §DD-11` and
`shell_vs_interior_diagnostic.md` report 6158 at K=3 (30,000 NFE). The 1362 figure is from
an early exploratory run. Cite run configuration explicitly; the DD-11 table is the reference.
`code_alignment_backlog.md Item 7` flags the CLAUDE.md update.

**FLAG B — Section outline version drift.**
`style_guide.md §3.1` lists a nine-section outline that still includes multi-site DRB sections.
Authoritative outline is in `manuscript_main_draft.md`. style_guide.md §3.1 is descriptive, not
prescriptive.

**FLAG C — "Hazard-space is emergent from parameter-space" nuance.**
MOEA-FIND searches parameter space to achieve hazard coverage — not the reverse. When migrating
any notes text into the manuscript, add this distinction. See `framing_anchor.md §4` and
`manuscript_main_draft.md §4`.

---

## literature/notes/ (quick-reference summaries)

`borgomeo_2015.md`, `parametric_generation.md`, `wheeler_2025.md`, `zaniolo_2024_FIND.md`

---

## scratch/ and archive/

**scratch/** — intro drafts, framing experiments, figure1 spec. Ephemeral; safe to delete.

**archive/** — `method_proposal.md` (KEEP-HISTORICAL; §3.1 architecture diagram and §8 parametric
reference list), `autonomous_plan.md`, `worktree_sync_audit.md`.
