# manuscript/ — Organization

*Established 2026-04-15. Authoritative on folder layout, file roles, and sweep policy.*

---

## Folder layout

```
manuscript/
├── *.md                  master docs — authoritative, continuously updated
├── ORGANIZATION.md       this file
├── INDEX.md              short pointer index to all master docs
├── figures/              figure scripts and generated outputs (outputs gitignored)
├── literature/
│   ├── *.pdf             PDFs (gitignored)
│   └── notes/            short summaries for quick drafting reference
├── scratch/              iterative and ephemeral work; safe to delete
└── archive/              KEEP-HISTORICAL files; do not edit
```

---

## Master docs (manuscript/ root)

One authoritative file per topic. Edit in place; do not fork into new files.

| File | Topic |
|---|---|
| `manuscript_main_draft.md` | Full manuscript prose — never delete |
| `supporting_info_draft.md` | Supplemental information draft — never delete |
| `framing_anchor.md` | DMDU / Herman Axis II / Moallemi Fork 3.2.2 framing reference |
| `reviewer_defenses.md` | Anticipated reviewer critiques and defensible responses |
| `publication_plan.md` | Phased work breakdown (Phases A–E) and exit criteria |
| `design_decisions.md` | Numbered design decisions DD-01 through DD-N |
| `code_alignment_backlog.md` | Code-to-manuscript alignment task list |
| `constraints_spec.md` | Constraint design, gap analysis, and calibration recommendations |
| `critical_analysis.md` | Honest risk assessment across eight topics |
| `research_questions.md` | RQ1–RQ8, figure sequence, SI plan |
| `shell_vs_interior_diagnostic.md` | K-ball interior coverage diagnostic results |
| `literature_review.md` | Structured literature review across five domains |
| `experiment_plan.md` | HPC deployment tracker and Sections 3.2 / 3.3 experiment design |
| `style_guide.md` | Manuscript prose rules — highest priority, overrides all |

---

## Other locations

**`literature/`** — PDFs stay gitignored where they are. `literature/notes/` holds short tracked
summaries (one per key paper) for quick lookup while drafting.

**`scratch/`** — Iterative work: intro drafts, framing experiments, figure specs under active
development. No retention guarantee. Content is safe to delete once merged into a master doc.
Do not cite scratch/ files as authoritative.

**`archive/`** — KEEP-HISTORICAL files marked with a banner. Do not edit. Retained for provenance
only.

**`figures/`** — Figure scripts (tracked) and generated outputs (gitignored).

---

## Policies

**Adding a new document.** First ask: is this a new authoritative topic, or is it exploratory?
- Authoritative → new root-level .md; add a row to INDEX.md.
- Exploratory → put it in scratch/.

**Promoting scratch → master.** Extract the load-bearing content and merge it into the relevant
master doc, then delete the scratch file.

**Archiving.** A file goes to archive/ when it is superseded but must be preserved for provenance.
Add a KEEP-HISTORICAL banner at the top before moving.

**Sweep policy.** At the start of any session, scratch/ files whose content has been absorbed into
a master doc may be deleted outright without ceremony.
