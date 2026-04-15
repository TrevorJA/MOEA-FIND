# Worktree Sync Audit
**Generated:** 2026-04-15  
**Auditor:** Claude Code (busy-shirley session)  
**Scope:** main + 6 overnight worktrees  
**Action taken:** read-only — no merges, no commits, no pushes

---

## 0. Executive Summary

**Critical finding:** No worktree branch has any commits ahead of `main`. All 6 branches sit at `5639f62` — identical to `main`. There is nothing to merge via git. The overnight sessions wrote to working-tree files (manuscript/*, scripts/*, src/*) that were never staged or committed.

Two classes of risk:

| Risk class | Where | Files at risk |
|---|---|---|
| **Uncommitted git-tracked files** | main checkout (unstaged mods + untracked new scripts) | 8 files, ~980 net lines |
| **Gitignored manuscript/ content** | zen-swanson (6 intro files), eager-dirac (2 files), main (autonomous_plan.md diverged) | 9 files |

`manuscript/` is in `.gitignore`. Everything under it is invisible to git — it cannot be merged, diffed via git, or recovered after a worktree prune.

---

## 1. Per-Worktree State

### main (`main`, commit `5639f62`)

**Git status — unstaged modifications (tracked files):**
```
M  scripts/10_plot_manuscript_figures.py      +338/-161  (substantial rewrite)
M  src/objectives.py                           +133/-0
M  src/plotting/analytic.py                   +433/-0    (near-rewrite)
M  src/plotting/trace_diagnostics.py          +237/-0
```
**Git status — untracked new files (not in .gitignore):**
```
?? scripts/diag_shell_vs_interior.py          (shell-vs-interior coverage diagnostic)
?? scripts/diag_shell_vs_interior_summary.py  (dimension-sweep summary figure)
?? src/plotting/architecture.py               (algorithmic schematics for manuscript)
?? src/plotting/map.py                        (DRB study-area map, Figure 4)
```
**Commits ahead of main:** 0 (this IS main)  
**Commits behind main:** n/a  
**Last 24 h commits:** none  
**Manuscript/ content:** ~40 files (all gitignored); includes `manuscript_main_draft.md`, `methods_audit.md`, `generator_constraints_review.md`, etc. These appear to be the primary destination where eloquent-northcutt and elastic-merkle sessions wrote their manuscript deliverables.

**Risk: HIGH.** ~980 lines of source-code work and 4 new scripts are not committed to any branch.

---

### zen-swanson — Intro Revision (`claude/zen-swanson`, `5639f62`)

**Git status:** clean (0 commits ahead, 0 behind)  
**Git-tracked changes:** none  
**Manuscript/ files (gitignored — only exist in this worktree):**
```
manuscript/intro_critique_round1.md
manuscript/intro_critique_round2.md
manuscript/introduction_framing_moves.md
manuscript/introduction_revised.md
manuscript/introduction_revised_round1.md
manuscript/introduction_revised_round2.md
```
**Collision with main:** None of these filenames exist in main's manuscript/. Safe to copy.  
**Risk: HIGH.** 6 intro-revision files exist ONLY in this worktree's working directory. They will be lost if the worktree is pruned.

---

### eager-dirac — HPC Deployment Prep (`claude/eager-dirac`, `5639f62`)

**Git status — untracked new files (not in .gitignore):**
```
?? configs/section32_cannonsville.yaml
?? configs/section33_pywr_drb.yaml
?? environment_hpc.yml
?? scripts/hpc/preflight_test.py
?? scripts/hpc/section32_cannonsville_library.slurm
?? scripts/hpc/section32_cannonsville_moea.py
?? scripts/hpc/section32_cannonsville_moea.slurm
?? scripts/hpc/section33_pywr_drb_batch.py
?? scripts/hpc/section33_pywr_drb_batch.slurm
```
**Manuscript/ files (gitignored):**
```
manuscript/autonomous_plan.md    ← DIVERGED from main (see §3)
manuscript/hpc_deployment_status.md  ← new, not in main
```
**Risk: MEDIUM.** 9 HPC files not committed. The `autonomous_plan.md` conflict is the only true content collision in this audit.

---

### adoring-elion — Original Intro Lit Review (`claude/adoring-elion`, `5639f62`)

**Git status:** clean  
**Manuscript/ files:** none  
**Risk: NONE.** Session produced no working-tree output. Either the session didn't execute or its outputs were written to the main checkout's manuscript/ directory.

---

### eloquent-northcutt — Methods / Section 3.1 (`claude/eloquent-northcutt`, `5639f62`)

**Git status:** clean  
**Manuscript/ files:** none  
**Risk: NONE.** Session deliverables (Generator description revision, MOEA setup expansion, D* caveat, `methods_audit.md`) appear in main's `manuscript/` per `autonomous_plan.md` — they were written to the main checkout, not this worktree.

---

### adoring-moser — Non-Drought Constraints Experiment (`claude/adoring-moser`, `5639f62`)

**Git status:** clean  
**Manuscript/ files:** none  
**Risk: NONE.** Session produced no working-tree output.

---

### elastic-merkle — Kirsch Constraints Review (`claude/elastic-merkle`, `5639f62`)

**Git status:** clean  
**Manuscript/ files:** none  
**Risk: NONE.** Consistent with "may have landed on main already" — `generator_constraints_review.md` is present in main's manuscript/. Confirm that file was elastic-merkle's deliverable before pruning this worktree.

---

## 2. Files Touched per Worktree (vs `main`)

Since no worktree has commits ahead of main, `git diff --name-only main...HEAD` returns empty for all. The table below reflects working-tree-level changes only (git-tracked unstaged + untracked).

| Worktree | scripts/ | src/ | manuscript/ (gitignored) |
|---|---|---|---|
| **main** | `10_plot_manuscript_figures.py` (M), `diag_shell_vs_interior.py` (U), `diag_shell_vs_interior_summary.py` (U) | `objectives.py` (M), `plotting/analytic.py` (M), `plotting/trace_diagnostics.py` (M), `plotting/architecture.py` (U), `plotting/map.py` (U) | see §1 |
| **zen-swanson** | — | — | 6 intro-revision files (U) |
| **eager-dirac** | `hpc/` 6 files (U) | — | `autonomous_plan.md` (diverged), `hpc_deployment_status.md` (U) |
| **adoring-elion** | — | — | — |
| **eloquent-northcutt** | — | — | — |
| **adoring-moser** | — | — | — |
| **elastic-merkle** | — | — | — |

M = modified tracked file; U = untracked new file

---

## 3. Collision Map

### Git-tracked collisions
**None.** No worktree has committed any file that conflicts with main.

### Manuscript/ collisions (gitignored — manual reconciliation required)

| File | In main | In worktree | Verdict |
|---|---|---|---|
| `manuscript/autonomous_plan.md` | Yes — intro revision checked ✓, HPC deployment unchecked | Yes (eager-dirac) — HPC deployment checked ✓, intro revision unchecked | **COLLISION.** Different tasks marked complete. Merge manually: both sets of checkmarks should be ✓. |
| `manuscript/hpc_deployment_status.md` | No | Yes (eager-dirac) | Copy eager-dirac → main |
| `manuscript/intro_critique_round1.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |
| `manuscript/intro_critique_round2.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |
| `manuscript/introduction_framing_moves.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |
| `manuscript/introduction_revised.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |
| `manuscript/introduction_revised_round1.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |
| `manuscript/introduction_revised_round2.md` | No | Yes (zen-swanson) | Copy zen-swanson → main |

---

## 4. Consolidation Recommendations

### Recommended action order

**Step 1 — Commit main checkout's source changes (no conflicts possible)**
```bash
cd /path/to/MOEA-FIND
git add src/objectives.py src/plotting/analytic.py src/plotting/trace_diagnostics.py \
        scripts/10_plot_manuscript_figures.py \
        scripts/diag_shell_vs_interior.py scripts/diag_shell_vs_interior_summary.py \
        src/plotting/architecture.py src/plotting/map.py
git commit -m "..."
```
These files touch no other worktree. Zero collision risk.

**Step 2 — Commit eager-dirac's HPC files (no conflicts with main)**
```bash
cd .claude/worktrees/eager-dirac
git add configs/ scripts/hpc/ environment_hpc.yml
git commit -m "..."
```
Then merge or cherry-pick to main.

**Step 3 — Copy zen-swanson's intro files to main manuscript/**
```bash
cp .claude/worktrees/zen-swanson/manuscript/introduction_revised*.md manuscript/
cp .claude/worktrees/zen-swanson/manuscript/intro_critique*.md manuscript/
cp .claude/worktrees/zen-swanson/manuscript/introduction_framing_moves.md manuscript/
```
No conflict — none of these filenames exist in main.

**Step 4 — Manually reconcile autonomous_plan.md**
Both main and eager-dirac modified `manuscript/autonomous_plan.md`. Merge the two checkmark states:
- Keep intro revision ✓ from main's copy
- Keep HPC deployment ✓ from eager-dirac's copy

**Step 5 — Copy eager-dirac's hpc_deployment_status.md to main**
```bash
cp .claude/worktrees/eager-dirac/manuscript/hpc_deployment_status.md manuscript/
```

**Step 6 — Prune empty worktrees (safe)**
adoring-elion, eloquent-northcutt, adoring-moser are completely empty. After verifying elastic-merkle's deliverable is confirmed in main's `generator_constraints_review.md`, these four can be pruned:
```bash
git worktree remove .claude/worktrees/adoring-elion
git worktree remove .claude/worktrees/eloquent-northcutt
git worktree remove .claude/worktrees/adoring-moser
git worktree remove .claude/worktrees/elastic-merkle
```
zen-swanson and eager-dirac: prune ONLY after Steps 3–5 are complete.

---

## 5. Orphaned / Discardable Changes

| Item | Location | Assessment |
|---|---|---|
| adoring-elion working tree | `.claude/worktrees/adoring-elion/` | Empty. Safe to prune. |
| eloquent-northcutt working tree | `.claude/worktrees/eloquent-northcutt/` | Empty. Session deliverables confirmed in main's manuscript/. Safe to prune. |
| adoring-moser working tree | `.claude/worktrees/adoring-moser/` | Empty. Safe to prune. |
| elastic-merkle working tree | `.claude/worktrees/elastic-merkle/` | Empty. Confirm `generator_constraints_review.md` is the session output, then prune. |

---

## 6. Main Branch Activity (last 24 h)

```
(no commits in last 24 hours)
```
Last commit: `5639f62 2026-04-13 Updates`  
Nothing has been silently committed.

---

*End of audit. No files were modified, no commits created, no merges performed.*
