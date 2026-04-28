# Supporting Information: MOEA-FIND

*SCAFFOLD DRAFT. Section headers only. Each subsection body is a placeholder
with a pointer to the underlying authoritative source. Prose will be
drafted only after the corresponding main-text section has been written
against final empirical results, so that the SI is written to support
specific main-text claims rather than anticipating them.*

*Companion to the main text. Each section is cross-referenced from a
specific location in the main text. Target length is approximately three
thousand words plus figures and tables.*

*Last updated: 2026-04-28 (workflow reorganization; stage→subsection mapping added).*

## Workflow stage to SI subsection mapping

Each method-consideration driver under [workflows/](../../workflows/) earns
SI text. The mapping below is the contract: every kept driver in stages 02
and 04, plus the dimension sweep in stage 01 and the verification/discovery
drivers in stages 06–07, is referenced below. Drivers that produce only
intermediate calibration data (no SI figures) are still tracked, with the
subsection that consumes them noted.

| Driver | SI subsection |
|---|---|
| `01_analytic_validation/eps_nfe_sweep.py`            | SI-3 (epsilon and NFE sensitivity) |
| `01_analytic_validation/dimension_sweep.py`          | SI-2 (interior-filling coverage) |
| `02_calibration/constraint_calibration.py`           | SI-4b (constraint regime ablation, hydrologic arm) |
| `02_calibration/dv_uniformity_calibration.py`        | SI-4b (constraint regime ablation, DV-space arm) |
| `02_calibration/wrapper_fidelity.py`                 | SI-5 (DV parameterisation, fidelity check) |
| `02_calibration/wrapper_geometry.py`                 | SI-5 (DV parameterisation, geometry check) |
| `02_calibration/kirsch_convergence.py`               | SI-4 (Borg convergence + Kirsch wall-clock) |
| `02_calibration/metric_blocks.py`                    | SI-6 (metric set + SSI sensitivity) |
| `04_moea_find_single_site/wrapper_mode_ablation.py`  | SI-5 (DV parameterisation ablation) |
| `04_moea_find_single_site/wrapper_mode_compare.py`   | SI-5 (DV parameterisation ablation) |
| `04_moea_find_single_site/dv_uniformity_ablation.py` | SI-4b (constraint regime ablation) |
| `04_moea_find_single_site/dv_uniformity_compare.py`  | SI-4b (constraint regime ablation) |
| `04_moea_find_single_site/event_level.py`            | SI-9 (event-level Kirsch — placeholder below) |
| `06_pywrdrb_reeval/verify_drought_coverage.py`       | SI-10 (Pareto coverage verification — placeholder below) |
| `07_scenario_discovery/satisficing_sweep.py`         | SI-11 (satisficing manifold + GBT — placeholder below) |
| `07_scenario_discovery/scenario_discovery_plots.py`  | SI-11 (satisficing manifold + GBT) |

Each subsection targets approximately 2-3 multi-panel figures (clean,
multi-panel academic style) plus supporting tables. Subsections SI-9,
SI-10, and SI-11 below are new placeholders added with this mapping.

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

> *Placeholder. Driver: `workflows/01_analytic_validation/eps_nfe_sweep.py`. SLURM
> array configuration in `workflows/0N_<stage>/slurm/`. Recommended default epsilon
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
> live in `workflows/0N_<stage>/` and `workflows/0N_<stage>/`; see
> `planning/code_state.md`.*

---

## SI-5. Sensitivity to the decision-variable parameterisation

*Cross-referenced from main-text Section 2.2.4.*

> *Placeholder. The decision-variable parameterisation of the Borg–generator
> coupling is not yet final. The production choice and any wrapper-mode
> ablation will be drafted here once settled. Current drivers live in
> `workflows/0N_<stage>/`; see `planning/code_state.md`.*

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

## SI-9. Event-level Kirsch objective formulation

*Cross-referenced from main-text Section 2.2 (Fig 5 inset).*

> *Placeholder. Event-level objective formulation that ties Kirsch DVs
> directly to drought-event characteristics (rather than aggregate
> series-level statistics). Driver:
> [`workflows/04_moea_find_single_site/event_level.py`](../../workflows/04_moea_find_single_site/event_level.py).
> Target 2-3 multi-panel figures: per-event Pareto, sensitivity to event
> definition window, comparison to series-level objectives. Pending HPC
> production runs.*

---

## SI-10. Pareto-archive drought coverage verification

*Cross-referenced from main-text Section 3 (Pywr-DRB re-evaluation
inputs).*

> *Placeholder. Verifies that the Pareto archive supplied to Pywr-DRB
> covers the drought-characteristic space without gaps before the
> expensive simulation is committed. Driver:
> [`workflows/06_pywrdrb_reeval/verify_drought_coverage.py`](../../workflows/06_pywrdrb_reeval/verify_drought_coverage.py).
> Target 2-3 multi-panel figures: drought-space coverage scatter,
> per-axis marginal histograms, FDC subset diagnostics.*

---

## SI-11. Satisficing manifold and GBT classifier diagnostics

*Cross-referenced from main-text Section 3.3 / Fig 9.*

> *Placeholder. Per-definition GBT decision-boundary overlays in
> drought-feature space, classifier ROC summary, manifest status table.
> Drivers:
> [`workflows/07_scenario_discovery/satisficing_sweep.py`](../../workflows/07_scenario_discovery/satisficing_sweep.py)
> and
> [`workflows/07_scenario_discovery/scenario_discovery_plots.py`](../../workflows/07_scenario_discovery/scenario_discovery_plots.py).
> Manifest of binary satisficing rules:
> [`workflows/07_scenario_discovery/satisficing_manifest.yaml`](../../workflows/07_scenario_discovery/satisficing_manifest.yaml).
> Target 2-3 multi-panel figures.*

---

> *Figure and table numbering in this Supporting Information will be
> finalised once the main-text figures are locked. Any SI section whose
> upstream output has not been produced by the submission date will be cut
> or replaced with a data-tables-only entry.*
