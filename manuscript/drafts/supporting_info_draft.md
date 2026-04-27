# Supporting Information: MOEA-FIND

*SCAFFOLD DRAFT. Section headers only. Each subsection body is a placeholder
with a pointer to the underlying authoritative source. Prose will be
drafted only after the corresponding main-text section has been written
against final empirical results, so that the SI is written to support
specific main-text claims rather than anticipating them.*

*Companion to the main text. Each section is cross-referenced from a
specific location in the main text. Target length is approximately three
thousand words plus figures and tables.*

*Last updated: 2026-04-27 (scaffold reset).*

---

## SI-1. Manhattan-distance auxiliary objective and the K-dimensional coverage argument

*Cross-referenced from main-text Section 2.2.3.*

> *Placeholder. The L1 formulation is fixed: $f_j = D_j$ for $j = 1, \ldots,
> K$ and $f_{K+1} = \lVert D - D^* \rVert_1$. The non-dominance argument and
> the codimension-one affine-subspace identity are documented in
> `governance/design_decisions.md` §DD-11 and implemented in
> `src/objectives.py`. The formal derivation will be drafted into this
> section once the empirical sections that motivate the framing exist.*

---

## SI-2. Empirical verification of interior-filling coverage

*Cross-referenced from main-text Section 3.1.*

> *Placeholder. Preliminary dimension-sweep evidence on a constrained
> analytic benchmark is documented in
> `evidence/shell_vs_interior_diagnostic.md` (EpsNSGAII stand-in). The
> benchmark, dimensionalities tested, and final tabulation are not yet
> settled and will be finalised once §2.3 of the main text is written and
> production Borg results on HPC are in hand.*

---

## SI-3. Epsilon and function-evaluation sensitivity on the analytic test problem

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. Driver: `workflows/experiments/03_eps_nfe_sweep.py`. SLURM
> array configuration in `workflows/slurm/`. Recommended default epsilon
> vector pending HPC sweep completion.*

---

## SI-4. Convergence diagnostics for Borg MOEA

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. Multi-seed hypervolume curves, archive size traces, and
> adaptive operator selection probabilities for every main-text experiment.
> Pending HPC production runs.*

---

## SI-4b. Constraint regime ablation

*Cross-referenced from main-text Section 2.2.4 and `governance/design_decisions.md` §DD-14.*

> *Placeholder. Constraint regime ablation. The set of regimes compared and
> the production choice are not yet final. Current drivers and diagnostics
> live in `workflows/experiments/` and `workflows/diagnostics/`; see
> `planning/code_state.md`.*

---

## SI-5. Sensitivity to the decision-variable parameterisation

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. The decision-variable parameterisation of the Borg–generator
> coupling is not yet final. The production choice and any wrapper-mode
> ablation will be drafted here once settled. Current drivers live in
> `workflows/experiments/`; see `planning/code_state.md`.*

---

## SI-6. Sensitivity to the drought metric set and SSI accumulation period

*Cross-referenced from main-text Section 2.1.*

> *Placeholder. The production metric set (DD-04) is the `primary`
> preset; sensitivity to the alternative presets (`extreme_event`,
> `trace_fdc`) and to the SSI accumulation period will be reported here
> once HPC results land. The metric infrastructure
> (`src/drought_metrics.py`) makes swapping metric sets a single CLI
> flag, so the ablations are inexpensive in code-engineering time even
> if HPC compute is gated.*

---

## SI-7. Per-site plausibility diagnostics (multi-site extension)

*Cross-referenced from main-text Section 2.4. Included only if the
multi-site extension makes the main text.*

> *Placeholder. Multi-site DRB extension is currently deferred to a
> follow-up paper per `planning/publication_plan.md` §4. This section will
> most likely be cut from the single-site submission. Site list and
> diagnostics are not finalised.*

---

## SI-8. Pareto archive reference tables

*Cross-referenced from main-text Section 3.2.*

> *Placeholder. Reference tables for the production archive will be
> released as supplementary files at submission. The exact columns released
> will be defined once §3.2 is written.*

---

> *Figure and table numbering in this Supporting Information will be
> finalised once the main-text figures are locked. Any SI section whose
> upstream output has not been produced by the submission date will be cut
> or replaced with a data-tables-only entry.*
