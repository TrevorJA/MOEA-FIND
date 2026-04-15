# Figure 1 Design Specification: Complementary Frameworks for Synthetic Scenario Design

*Written 2026-04-14. Purpose: specify the panel-by-panel layout, visual grammar,
citation anchors, and design rationale for the manuscript's opening conceptual
figure. This spec is grounded in a literature survey of schematic conventions used
in the DMDU, exploratory modeling, and synthetic hydrology communities. No figure
image is generated here — this is a design document for the illustrator/coding
phase.*

---

## 1. Literature survey: schematic conventions found

The table below records which papers contain a schematic figure showing a data-flow
or workflow across two or more phases of analysis, what spaces those figures show,
and the visual grammar each paper uses. Sources were read via PDF extraction and
web search on 2026-04-14.

| Paper | Schematic figure? | Spaces shown | Visual grammar | Notes |
|---|---|---|---|---|
| Bonham et al. (2024) | Yes — Fig. 1 (confirmed caption) | Uncertainty space (scenarios), decision space (policies), outcome space (performance + rankings) | Three labeled rectangular boxes connected by downward arrows; each box corresponds to one step; step labels in bold above each box | Caption: "Framework overview. In step 1... step 2... step 3..." — purely step-sequential, no space annotations |
| Bonham et al. (2025) | Yes — Fig. 1 (confirmed caption) | SOW ensemble space, simulation output space, outcome/label space, uncertain-factor space | Four boxes in a vertical or L-shaped flow: SOWs+Policies → Simulation → Performance metrics → Decision-relevant outcomes → Factor mapping | The boxes are explicitly labeled with the variable type (SOW, metric, outcome, factor); this is the cleanest published template for the DMDU pipeline |
| Hadjimichael et al. (2020, Earth's Future) | Yes — schematic described in text; explicit panel I–III structure | Parameter/uncertainty space (14-dim LHS), simulation output space, performance/robustness space | Three panel structure: panel I (factor definitions), panel II (performance across realizations), panel III (user performance across scenarios) | The three-panel structure maps directly to Moallemi's three spaces |
| Moallemi et al. (2020) | Implied — encoding error prevented full extraction; abstract describes four sequential phases | Decision space, uncertainty space, outcome space (explicitly named as "three spaces") | Four-phase horizontal flow (problem context → decision framing → solution evaluation → recommendations); "exploratory modeling" populates the uncertainty-to-outcome arrow | THE vocabulary anchor for the three-space taxonomy; the three-space framing is the direct citation basis for Figure 1's shared axis |
| Brown et al. (2012) | Yes — Fig. 1 (confirmed, three-step process schematic) | Climate forcing space (T, P grid), performance/outcome space (climate response function), probability space (climate projections) | Three labeled steps as a horizontal flow with a threshold boundary drawn across the climate space; step 3 overlays climate projections on the boundary | The "climate response function" is the key visual object: a grid of climate states colored by performance outcome, with a failure boundary contour. This is the direct ancestor of MOEA-FIND's feasible-region concept. |
| Herman et al. (2015) | Table 1 (taxonomy table, not a flowchart) | Four methodological dimensions: alternative generation, SOW sampling, robustness quantification, sensitivity analysis | A comparison table across four frameworks (RDM, Decision Scaling, Info-Gap, MORDM) — not a data-flow schematic | Still useful: establishes the four dimensions that differentiate frameworks and could inform the annotation layer of Figure 1 |
| Culley et al. (2016) | Yes (implied by scenario-neutral methodology) | Climate exposure space (forcing perturbations), optimization space (policy search), adaptive capacity boundary | Shows an "exposure space" with a boundary curve beyond which no operating rule can maintain performance | Closest conceptual ancestor to MOEA-FIND panel: defines the feasible region boundary and shows where it lies in the forcing space |
| Weaver et al. (2013) | Yes — Fig. 2 (confirmed as "Steps in a robust decision-making analysis," reprinted from Lempert and Groves 2010) | XLRM space (uncertainties X, levers L, relationships R, measures M), SOW ensemble, strategy space, vulnerability region | Five-step sequential flowchart: problem framing → candidate strategies → ensemble evaluation → scenario discovery (PRIM) → tradeoff analysis | The canonical RDM workflow figure; used as Figure 2 in multiple subsequent papers; widely recognized in the DMDU community |
| Maier et al. (2016) | Yes — Fig. 1 (confirmed as three modeling paradigms) | Three paradigms: predict-then-act (single trajectory), probabilistic ensemble, exploratory scenario (deeply uncertain) | Three parallel panels or regions showing distinct epistemic stances toward the future; transitions between paradigms shown as arrows | Framework-comparison figure rather than data-flow; useful for positioning argument but not a template for the data-flow schematic |
| Marchau et al. (2019) (Springer handbook, open access) | Yes — multiple workflow figures, one per method chapter; overview comparison in Ch. 1 | XLRM for RDM; adaptation tipping points for DAPP; decision branches for DAP; possibility space for Info-Gap | Each method has a dedicated linear or branching flowchart; Ch. 1 overview figure (Fig. 1.3) places all methods in a comparative matrix | The most thorough source for DMDU workflow figures; RDM workflow is the best-documented and most widely cited |

**Key observations from the survey:**

1. The dominant convention is a **left-to-right or top-to-bottom sequential flow** of
   rectangular or rounded boxes, one per phase. The phases are: (sample/perturb
   inputs) → (simulate model) → (evaluate performance) → (analyze/discover).

2. All frameworks operate in the **three spaces** described by Moallemi et al.
   (2020): uncertainty/forcing space, simulation space, and outcome/performance
   space. The primary design choice in each framework is *where* the structured
   sampling effort is concentrated.

3. The **Brown (2012) decision scaling figure** is the most visually distinctive
   because it draws the threshold boundary as a contour in the forcing space,
   making the feasible-region concept visible. This is the convention that maps
   most directly to MOEA-FIND.

4. The **Bonham (2025) Fig. 1** is the best published template for the DMDU
   pipeline: four labeled boxes (SOWs + policies → simulation → metrics →
   outcomes → factor mapping) with explicit variable-type labels. The MOEA-FIND
   Figure 1 should use the same box-and-arrow grammar adapted to four frameworks
   rather than four phases.

5. No published paper places four frameworks side by side in a single figure using
   consistent visual grammar across all four panels. The closest precedent is
   Herman et al. (2015) Table 1, which uses a comparison table rather than
   flowcharts. **Figure 1 fills a gap in the literature's own self-representation.**

---

## 2. Design rationale

### 2.1 Main message

MOEA-FIND is a complementary fourth approach to synthetic scenario design, not a
replacement for existing methods. The figure makes this visible by showing all four
frameworks with the same boxes and arrows but with the structured design step
(the "gold star" step) applied at different points in the flow. The reader can
see in one glance that traditional methods design for coverage upstream of the
simulation model (in the forcing variable or parameter space), while MOEA-FIND
designs for coverage downstream (in the outcome/hazard space), and that both are
valid depending on the research question.

### 2.2 The reversed-arrow device

The single most important visual element is the **reversed search arrow** in the
MOEA-FIND panel. In panels A, B, and C, information flows strictly
left-to-right: inputs → generator → scenarios → model → outcomes. In panel D,
an additional arrow runs from the desired-coverage specification in outcome space
backward through the optimizer to the generator inputs, and the optimizer's
output drives the generator. This reversed flow is the visual manifestation of
"hazard-space targeting" — the reader's eye is disrupted at exactly the right
conceptual moment.

### 2.3 The shared hazard-space strip

A horizontal strip spanning all four panels at the bottom of the figure shows
the qualitative coverage pattern that each framework produces in drought hazard
space, using the same axis system (e.g., mean severity on the horizontal axis,
mean duration on the vertical axis). This shared axis makes explicit that all
four approaches are ultimately evaluated against the same hazard space, but only
MOEA-FIND controls coverage in that space by construction. The strip requires no
caption text to communicate this — the reader sees it directly.

### 2.4 Positioning MOEA-FIND last

Placing MOEA-FIND as the rightmost panel (panel D) exploits the reader's
left-to-right scanning direction. Panels A, B, and C build the expectation of
a forward flow; panel D breaks it with the reversed arrow. This is more
communicative than placing MOEA-FIND first or second, where the deviation from
convention would not yet have a baseline to violate.

### 2.5 Consistent grammar prevents false novelty claims

Using the same boxes, shapes, and arrow types across all four panels
communicates that MOEA-FIND uses existing components (historical data, a
streamflow generator, a simulation model, robustness analysis) assembled in a
different order. This reinforces the "complementary" framing and pre-empts the
"this is just a rebranding of multi-objective optimization" critique by showing
that the novelty is specifically in the *direction* of the structured design step,
not in the components themselves.

---

## 3. Panel-by-panel layout specification

**Overall dimensions:** Four panels arranged in a single horizontal row, each
approximately 40 mm wide × 90 mm tall, with 5 mm gaps between panels.
Total figure width: 175 mm (single-column WRR width). Total height: 90 mm.

The shared hazard-space strip occupies an additional 20 mm below the four panels,
spanning the full 175 mm width.

**Shared box-and-arrow grammar (all panels):**

| Element | Shape | Color (hex) | Meaning |
|---|---|---|---|
| Forcing/input space | Rounded rectangle | #4472C4 (steel blue) | Historical data, parameter ranges, climate inputs |
| Streamflow/scenario generator | Rectangle, thin border | #ED7D31 (orange) | Kirsch bootstrap, climate perturbation operator, LHS sampler |
| Scenario ensemble | Stack of three overlapping rectangles | #A5A5A5 (medium gray) | The collection of synthetic traces or SOWs |
| Simulation model | Rectangle, thick border | #ED7D31 (orange) | Hydrological or water supply model (Pywr-DRB, StateMod) |
| Outcome/hazard space | Rounded rectangle | #548235 (forest green) | Performance metrics, drought characteristics |
| Analysis step | Rounded rectangle | #8064A2 (light purple) | Scenario discovery, sensitivity analysis, robustness |
| Structured design marker | Gold arrow label "DESIGN" | #FFC000 (gold) | Where the analyst controls coverage |
| MOEA optimizer box | Rectangle, thick red border | #C00000 (deep red) | MOEA-FIND specific; absent in other panels |
| Forward arrows | Solid black, 1 pt | #000000 | Information flow from left to right |
| Search/optimization arrow | Solid red, 1.5 pt, dashed | #C00000 | Backward flow from hazard space to generator (panel D only) |

**Shared hazard-space strip:** Same green (#548235) background, 20 mm tall,
spanning all four panels. Each panel's region shows a dot cloud representing the
projected coverage that framework produces in the (severity, duration) plane.
Dot positions are schematic (not from actual data) and drawn to communicate the
qualitative coverage pattern described for each panel below.

---

### Panel A — Forward exploratory modeling

**Citation anchor:** Hadjimichael et al. (2020), Moallemi et al. (2020),
Quinn et al. (2020), Bonham et al. (2025)

**Title text (above panel):** "Forward exploratory modeling"

**Box layout (top to bottom):**
1. Blue rounded rectangle: "Historical forcing data / Uncertain parameters"
   - Sub-label (small italic): "X: uncertainties"
2. Arrow down: labeled "LHS design" (gold "DESIGN" label on this arrow)
3. Orange rectangle: "Streamflow generator / model"
4. Arrow down
5. Gray stacked rectangles: "Large ensemble of SOWs"
6. Arrow down
7. Orange thick-border rectangle: "Simulation model (Pywr-DRB)"
8. Arrow down
9. Green rounded rectangle: "Performance metrics (reliability, shortage)"
10. Arrow down
11. Purple rounded rectangle: "Scenario discovery in parameter space (PRIM / GBT)"

**Small annotation in lower-left of panel:**
"Coverage of hazard space: emergent"

**Shared strip (panel A region):** Scattered dots distributed non-uniformly,
clustered toward moderate severity and moderate duration, sparse at the
high-severity high-duration corner. Annotation: "uncontrolled."

---

### Panel B — Bottom-up vulnerability analysis

**Citation anchor:** Brown et al. (2012), Herman et al. (2016), Culley et al. (2016)

**Title text:** "Bottom-up vulnerability / stress testing"

**Box layout (top to bottom):**
1. Blue rounded rectangle: "Historical forcing data / Climate perturbation range"
   - Sub-label: "Forcing variable grid (e.g., ΔP, ΔT)"
2. Arrow down: labeled "Systematic grid or LHS over forcing space" (gold "DESIGN" label)
3. Orange rectangle: "Weather generator / perturbed traces"
4. Arrow down
5. Gray stacked rectangles: "Perturbed scenario ensemble"
6. Arrow down
7. Orange thick-border rectangle: "Simulation model"
8. Arrow down
9. Green rounded rectangle: "Performance outcomes (climate response function)"
   - Inset: small box-within-box showing threshold boundary as a dashed contour
     dividing green box into "acceptable" and "failure" regions
10. Arrow down
11. Purple rounded rectangle: "Failure boundary in forcing space"

**Additional dashed feedback arrow:** From the green "Performance outcomes" box
back up to the "Perturbed scenario ensemble" box, labeled "identifies where
failures occur in forcing space." This is the Brown/Culley direction: start from
outcomes to refine input sampling near the threshold.

**Small annotation:** "Coverage of hazard space: threshold-focused"

**Shared strip (panel B region):** Dots concentrated near a diagonal threshold
line, sparse in the low-severity / low-duration corner and sparse in the extreme
corner. Annotation: "boundary-focused."

---

### Panel C — Library-and-subsample

**Citation anchor:** Bonham et al. (2024), Herman et al. (2016)

**Title text:** "Library-and-subsample"

**Box layout (top to bottom):**
1. Blue rounded rectangle: "Historical streamflow record"
2. Arrow down
3. Orange rectangle: "Kirsch-Nowak bootstrap generator (unsteered)"
4. Arrow down
5. Gray stacked rectangles: "Large pre-generated library (10,000+ traces)"
   - Sub-label: "covers input space broadly"
6. Arrow down: labeled "cLHS subsampling in input space" (gold "DESIGN" label)
7. Smaller gray stacked rectangles: "Efficient subsample (N traces)"
   - Sub-label: "covers input space uniformly"
8. Arrow down
9. Orange thick-border rectangle: "Simulation model"
10. Arrow down
11. Green rounded rectangle: "Performance metrics"
   - Sub-label: "hazard-space coverage: indirect"
12. Arrow down
13. Purple rounded rectangle: "Robustness ranking / scenario discovery"

**Small annotation:** "Coverage of hazard space: input-space projection"

**Shared strip (panel C region):** Dots somewhat more uniformly distributed than
panel A but still showing nonlinear projection distortion — moderate clustering in
center, thin tails at extremes. Annotation: "input-projected."

---

### Panel D — MOEA-FIND (this paper)

**Citation anchor:** This paper

**Title text:** "MOEA-FIND: hazard-space targeting"

**Box layout (top to bottom):**
1. Blue rounded rectangle: "Historical streamflow record"
2. Arrow DOWN from blue box to orange generator (same as other panels)
3. Orange rectangle: "Kirsch-Nowak bootstrap generator"
4. Arrow DOWN from generator: "produces candidate trace"
5. Green rounded rectangle: "Drought hazard characteristics (D₁, D₂, D₃)"
   - This box is prominent — same size as blue box
6. Arrow DOWN from green box to purple analysis
7. Purple rounded rectangle: "Feasible hazard region discovered; structured ensemble"

**NEW elements specific to panel D:**
- Deep red rectangle (MOEA optimizer) positioned between the green "Drought
  hazard" box and the blue "Historical forcing" box, to the LEFT of the main
  vertical flow: "Borg MOEA: Manhattan-distance objective + ε-dominance archive"
- **RED DASHED ARROW** from "Drought hazard characteristics" box LEFT and UP to
  the red MOEA box: labeled "desired coverage in hazard space" — this is the
  reversed arrow
- **RED SOLID ARROW** from the MOEA box DOWN to the orange generator box:
  labeled "steered decision variables" — MOEA drives the generator
- Gold "DESIGN" label placed on the red dashed arrow
- Small annotation BELOW the green hazard box: "anti-ideal placement → feasible
  region boundary"

**Panel D has no LHS step, no subsampling step, and no thick-border simulation
model box** — the simulation role is played by the generator + drought metric
extraction pipeline. The Pywr-DRB simulation used in Section 3.3 for vulnerability
demonstration is not part of the core MOEA-FIND loop and is NOT shown in this panel;
the panel shows only the ensemble generation step.

**Small annotation:** "Coverage of hazard space: structured by construction"

**Shared strip (panel D region):** Dots arranged in a near-uniform grid pattern
filling the feasible region, including the interior. Annotation: "structured."

---

## 4. Caption draft

*Following WRR house style: starts with the object of the figure, uses
"Panel (a)... Panel (b)..." enumeration, no bold leads, no results previewed.*

Figure 1. Schematic comparison of four complementary approaches to synthetic
scenario ensemble design for water supply vulnerability analysis. Each panel
follows the same visual grammar: blue boxes represent historical forcing data and
uncertain parameter ranges, orange boxes represent streamflow generators and
simulation models, green boxes represent performance outcomes and drought hazard
characteristics, and purple boxes represent analysis and scenario discovery steps.
Solid black arrows indicate forward information flow; the red dashed arrow in
panel (d) indicates the direction of the structured design effort. The gold
"DESIGN" label marks the step at which each approach controls the coverage
properties of the resulting ensemble. The horizontal strip at the bottom of each
panel shows the qualitative pattern of coverage produced in the drought hazard
space of interest (drought severity on the horizontal axis, drought duration on
the vertical axis). Panel (a): Forward exploratory modeling applies a Latin
hypercube design in the uncertain parameter space and propagates scenarios
forward through the simulation model; coverage of the drought hazard space is
an emergent property of the generator map (Hadjimichael et al., 2020; Moallemi
et al., 2020; Quinn et al., 2020). Panel (b): Bottom-up vulnerability analysis
applies a systematic grid or space-filling design in the forcing variable space
and maps the threshold boundary in that space using the simulated performance
outcomes (Brown et al., 2012; Herman et al., 2016; Culley et al., 2016). Panel
(c): Library-and-subsample applies a space-filling subsampling criterion in the
generator input coordinates of a pre-generated scenario library (Bonham et al.,
2024); coverage of the drought hazard space is an indirect consequence of the
input-space design. Panel (d): MOEA-FIND (this paper) reverses the design
direction: a multi-objective evolutionary algorithm steers the Kirsch-Nowak
generator specifically to produce synthetic traces whose drought characteristics
span the feasible drought hazard region, and the convergence of the
epsilon-dominance archive is the structured ensemble.

---

## 5. Inspiration sources and what was drawn from each

**Bonham et al. (2024), Fig. 1.** The three-step sequential box-and-arrow
layout with step labels above each box is the direct template for the panel
structure. The explicit labeling of what gets calculated at each step ("robustness
is calculated," "space-filling statistics are calculated") informs the annotation
style inside boxes. The performance objectives described in their Fig. 2 caption
(reliability, resiliency, vulnerability as green/orange/purple) directly inform
the metric definitions used in the Section 3.3 redesign.

**Bonham et al. (2025), Fig. 1.** The explicit naming of variable types
(SOWs, policies, performance metrics, decision-relevant outcomes, uncertain
factors) at each box informs the sub-label conventions in each box in
Figure 1 panels A through C. The four-box pipeline (SOW ensemble →
simulation model → performance metrics → decision-relevant outcomes →
factor mapping) is the direct template for panels A and C.

**Brown et al. (2012), Fig. 1.** The "climate response function" concept —
showing performance outcomes as a colored region in the forcing space with a
threshold boundary contour — is the template for the inset in panel B's
outcome box and for the shared hazard-space strip at the bottom of all panels.
The three-step structure of Decision Scaling (threshold identification →
stochastic analysis → risk estimation) also informs the three-component
bottom-up vulnerability flow in panel B.

**Weaver et al. (2013), Fig. 2.** The RDM workflow (problem framing → generate
strategies → stress-test across SOW ensemble → scenario discovery → tradeoff
analysis) is the template for the forward-flow arrow structure common to panels
A, B, and C. The explicit labeling of the scenario discovery step as distinct
from the simulation step informs the placement of the purple analysis box below
the green outcome box.

**Hadjimichael et al. (2020) three-panel structure.** Panel I (factor
definitions) / Panel II (performance across realizations) / Panel III (user
performance across scenarios) maps to the blue / orange-gray / green portions of
each panel in Figure 1, confirming that the three-space decomposition is
recognizable to the target audience.

**Visual grammar summary:** The figure adopts the box-and-arrow convention
universal in the sampled literature, with rounded rectangles for "spaces" and
sharp rectangles for "processes," consistent with Bonham (2024, 2025) and
Weaver (2013). The color coding (blue = inputs, orange = processes, green =
outcomes, purple = analysis) is new to this figure but maps to intuitive
associations (cool = known/controlled, warm = computed, green = desired
outcome, purple = insight). The reversed red dashed arrow is a deliberate
departure from all sampled predecessors; no prior paper in the survey uses
a backward arrow in this type of schematic, making it immediately visible
as the distinctive element of the MOEA-FIND panel.

---

## 6. Open design questions for author review

1. **Three vs. four panels.** The current spec has four panels (A: forward
   exploratory modeling, B: bottom-up vulnerability, C: library-and-subsample,
   D: MOEA-FIND). If four panels make the figure too narrow at 175 mm, a three-
   panel version could merge A and C into a single "input-space sampling" panel
   and treat B and D as the two contrasting alternatives. The four-panel version
   is preferred because the library-and-subsample (Bonham 2024) is the direct
   methodological comparator and deserves its own panel.

2. **Simulation model box in panel D.** The current spec omits the Pywr-DRB
   simulation model from panel D because it is not part of the MOEA-FIND loop.
   An alternative is to include it as a dashed-border box to the right of the
   green hazard box, labeled "downstream simulation (optional)" with an arrow
   from the ensemble to the model. This would make the figure show the
   complete workflow including the Section 3.3 vulnerability demonstration.
   Recommended: include as dashed-border optional box for completeness.

3. **Hazard-space strip dot patterns.** The qualitative dot patterns in the
   shared strip are schematic, not data-derived. Once the Cannonsville results
   are available, replace the schematic dots with actual projected scatter plots
   from the three baselines and the MOEA-FIND archive for the most compelling
   version of the figure.

4. **Title text for each panel.** The four title texts ("Forward exploratory
   modeling," "Bottom-up vulnerability / stress testing," "Library-and-subsample,"
   "MOEA-FIND: hazard-space targeting") are informal labels for design purposes.
   Final panel titles should follow the WRR house style for multi-panel figures
   (typically just "(a)", "(b)", etc. in the figure itself, with full descriptions
   in the caption).
