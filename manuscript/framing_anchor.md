# Framing Anchor Synthesis

**Purpose:** Reading-and-synthesis pass over Herman et al. (2015), Moallemi et al. (2020a, 2020b), and
four canonical robustness-framework figures. No Introduction text is drafted here. This document is a
reference to be consulted during the rewrite.

---

## 1. Herman (2015) Taxonomy Summary

**Full citation:** Herman, J.D., Reed, P.M., Zeff, H.B., and Characklis, G.W. (2015). How Should
Robustness Be Defined for Water Systems Planning under Change? *Journal of Water Resources Planning and
Management*, 141(10), 04015012.

### The Four Axes

Herman proposes a taxonomy of robustness frameworks structured along four independent axes. Exact language
from the paper's abstract and Figure 1 caption:

> "We propose a taxonomy of robustness frameworks to compare and contrast these approaches based on their
> methods of (1) alternative generation, (2) sampling of states of the world, (3) quantification of
> robustness measures, and (4) sensitivity analysis to identify important uncertainties."

**Axis I — Alternatives**
Methodological choices within each sub-category (from Fig. 1):
- Prespecified | Search
- Optimization | Design of experiments
- Single or Multi-Objective? | Well-characterized or deep uncertainty?

**Axis II — States of the World** (most relevant axis for MOEA-FIND)
Sub-choices listed in Fig. 1:
- Key Factors Assumed | Key Factors Discovered
- Prespecified Scenarios | Design of experiments
- Sampling outward from expected state, or global? | Well-characterized or deep uncertainty?

Herman's description of this axis:

> "Whether alternatives are prespecified or optimized, robustness frameworks typically proceed by
> evaluating their performance across a set of uncertain states of the world... One major aim of
> robustness analyses is to isolate the deeply uncertain factors responsible for system
> vulnerabilities. If decision makers have reason to believe a priori that a certain subset of
> uncertain factors will be most relevant, these may be prespecified. Alternatively, many uncertain
> factors may be sampled... Once the uncertain factors are identified, specific combinations of them
> can be sampled to create states of the world."

Key language: "the selection and sampling of the uncertain factors which make up plausible future states
of the world" (Fig. 1 caption, verbatim).

Herman describes two sampling modes for states of the world:
- *Global exploratory sampling* (Decision Scaling, RDM, MORDM): "an exploratory global sample of
  uncertain factors over plausible ranges"
- *Radial sampling from expected state* (Info-Gap): "samples radially outward from the expected future
  state in the uncertainty space until it encounters a state that causes failure"

**Axis III — Robustness Measures**
Sub-choices from Fig. 1:
- Expected Value | Satisficing | Regret
- Domain Criterion or Uncertainty Horizon? | Univariate or Multivariate Thresholds?
- Deviation from Best or Baseline?

Herman quotes Lempert and Collins (2007) for the overarching categories of regret and satisficing, and
Starr (1962) / Schneller and Sphicas (1983) for the *domain criterion*. Key passage:

> "The domain criterion quantifies the volume of the uncertain factor space in which a solution meets
> the decision makers' performance requirements."

**Axis IV — Robustness Controls**
Sub-choices from Fig. 1:
- Identify most sensitive factors | No sensitivity analysis
- Ranges | Ranking | Local OAT

Herman explicitly labels this as "consequence-oriented sensitivity analysis" encompassing *factor
mapping* (identifying sensitive ranges — PRIM) and *factor prioritization* (ranking all factors —
Sobol).

### Where MOEA-FIND Sits in the Herman Taxonomy

MOEA-FIND is a method operating **within Axis II (States of the World)** — specifically on the
*generation* of states of the world, not the generation of decision alternatives.

The methodological gap in Herman's taxonomy: Axis II lists "Design of experiments" (i.e., parameter
sampling) as the search-based option for generating states of the world. Herman does not contemplate
directed multi-objective search *in outcome/hazard space* as a scenario-generation mode. MOEA-FIND
fills this gap. It is a search-based approach to generating states of the world, but the search operates
in drought characteristic (hazard outcome) space rather than in the model parameter space. The sub-choice
that best approximates MOEA-FIND in Herman's language is: Search | Key Factors Discovered | Global
(coverage of the reachable hazard space) | Deep uncertainty. However, MOEA-FIND's search variable is
drought outcome, not parameter values — this is the taxonomic novelty.

**Most important exact phrase to adopt from Herman:** "the selection and sampling of the uncertain
factors which make up plausible future states of the world" — use this to introduce Axis II in the
Introduction, then argue that prior sampling approaches have operated in the *parameter* portion of this
space, whereas MOEA-FIND samples the emergent *hazard outcome* space.

---

## 2. Moallemi (2020) Workflow Vocabulary

### Paper A: Exploratory Modeling Taxonomy

**Full citation:** Moallemi, E.A., Kwakkel, J., de Haan, F.J., and Bryan, B.A. (2020). Exploratory
modeling for analyzing coupled human-natural systems under uncertainty. *Global Environmental Change*,
65, 102186.

**Core taxonomy:** Five approaches in Boxes 1–2. Each approach is characterized along three named
components:

| Component | Definition (Moallemi's exact language) |
|---|---|
| **Decision specification** | "how decision alternatives are generated" |
| **Uncertainty characterization** | "how scenarios are generated" |
| **Type of outcome implication** | "the inferences obtained from the analysis" |

**The five approaches and their uncertainty characterization modes:**

a) *Design of experiments* — "random values from the uncertainty space are systematically sampled to
generate/represent many possible scenarios/states of the world"

b) *Stress-testing* — "random values are systematically sampled from the uncertainty space to
generate/represent many possible scenarios. The generated scenarios are then investigated with subspace
partitioning... for stress-testing and scenario discovery"

c) *Worst-case scenario discovery* — "the uncertainty space is searched for specific scenarios that
could lead to the worst outcome"

d) *Many-objective optimization* — "a set of pre-specified scenarios can be selected a priori for the
evaluation of the decision alternatives" (uncertainty characterization is passive here — search is on the
decision side)

e) *Many-objective robust optimization* — "random values are systematically sampled from the
uncertainty space to generate/represent many possible scenarios for testing the robustness of the
decision alternatives"

**Critical observation:** All five of Moallemi's approaches characterize uncertainty through *parameter
sampling* (approaches a, b, d, e) or worst-case parameter search (c). None characterize uncertainty
through directed search in hazard outcome space with structured coverage objectives. MOEA-FIND is a
sixth mode not captured by this taxonomy — a directed, coverage-optimizing search within the
*outcome/hazard space* emergent from the parameter space.

**Moallemi's definition of a scenario** (verbatim, Global Environmental Change paper):

> "A scenario, in the context of exploratory modeling, is a fully specified realization of sampled
> values from the model parameter uncertainty space, representing a future state of the world."

This definition is anchored to *parameter sampling*. MOEA-FIND relaxes this: a scenario can be a
fully specified realization selected for structured coverage of *outcome/hazard* space, with the
parameter values as implicit inputs rather than the explicit sampling variable.

**Moallemi's description of the two modes (open exploration vs. directed search):**

> "Open exploration through the design of experiments and stress-testing is one way of investigating
> global properties of the assumption set... The weakness of open exploration emerges for systems with
> complex combinations of decisions that require the attainment of high levels of sustained performance
> often for conflicting objectives."

MOEA-FIND is a directed search method applied not to the decision space but to the uncertainty/scenario
space — a combination not explicitly discussed in Moallemi (2020a).

### Paper B: Structuring and Evaluating Decision Support Processes

**Full citation:** Moallemi, E.A., Zare, F., Reed, P.M., Elsawah, S., Ryan, M.J., and Bryan, B.A.
(2020). Structuring and evaluating decision support processes to enhance the robustness of complex
human–natural systems. *Environmental Modelling and Software*, 123, 104551.

**Four named stages (adapted from Tsoukias, 2008):**

1. **Representation of the problem context** — setting the stage: stakeholder analysis, clarifying
   current knowledge and assumptions, defining analytical goal(s)
2. **Problem framing** — formulating decisions, framing future scenarios, selecting robustness measures
3. **Evaluation and implementation** — developing or selecting an evaluation model, running
   computational experiments
4. **Decision recommendations and monitoring**

**Decision forks within each stage** (from Section 3, cited from Fig. 3 of the paper — 11 named forks):

Stage I forks: (3.1.1) Analysing stakeholders; (3.1.2) Clarifying current knowledge and assumptions;
(3.1.3) Defining analytical goal(s)

Stage II forks: (3.2.1) Formulating decisions; (3.2.2) Framing future scenarios; (3.2.3) Selecting
robustness measures

Stage III forks: (3.3.1) Developing or selecting an evaluation model; [additional forks to 3.3.2–3.3.3
not fully enumerated in read sections]

Stage IV forks: decision recommendations and monitoring

**Where MOEA-FIND sits (Moallemi structuring paper):** MOEA-FIND operates at Stage II, fork 3.2.2:
**Framing future scenarios**. From the paper:

> "There are two main aspects of framing future scenarios (i.e., states of the world) over which
> decisions are to be evaluated. The first is the identification of influential but uncertain factors
> shaping future scenarios... via participatory processes with stakeholders. The second is the means
> by which uncertain factors are selected for representation in the model-based quantitative assessment."

The paper identifies two modes for generating scenarios (states of the world):
- *Pre-specified standardised (and deterministic) frames of the future* (e.g., expert-defined scenarios)
- *Globally sampled stochastic analysis* (e.g., LHS, Monte Carlo over parameter uncertainty space)

MOEA-FIND introduces a third mode: directed search over hazard outcome space with structured coverage
objectives, where the scenarios are selected a posteriori based on their drought characteristics rather
than sampled uniformly from the parameter space. The paper does not enumerate this third mode.

**Taxonomy of robustness frameworks (Moallemi structuring, Section 2.2, citing Kwakkel and Haasnoot
2019 and Herman et al. 2015):**

> "Generation of decisions (i.e., alternatives): Candidate decisions can be pre-specified standardised
> alternatives... or can be generated through computational search..."
>
> "Generation of scenarios (i.e., states of the world): Scenarios can be specified as pre-specified
> standardised scenarios... or can be generated through exploratory modelling over the uncertainty space."
>
> "Measurement of performance (i.e., robustness): Different measures can be used to represent the
> performance of candidate decisions under uncertainty... Among them are descriptive statistical
> measures, satisficing measures, and regret measures."
>
> "Vulnerability analysis (i.e., robustness controls): Decision support can use different ways to
> isolate scenarios responsible for a specific behaviour of interest..."

**Key path-dependency insight:** Moallemi (2020b) emphasizes that choices at decision forks create
cumulative impacts ("path dependency"). The choice of scenario-generation method (fork 3.2.2) strongly
influences what robustness insights can be drawn. MOEA-FIND argues that parameter-space sampling (the
dominant mode) may not produce structurally complete coverage of the hazard outcome space — a
path-dependency concern about scenario generation.

---

## 3. Figure Grammar Inventory

### Figure 1: `herman_etal_2015_robustness_taxonomy.jpg`

**Source:** Herman et al. (2015), Figure 1.

**Visual grammar:**
- Four horizontally stacked rectangular panels, each shaded in a different pale tone (light blue-green
  for I and II; pale blue for III; white-gray for IV)
- Each panel has a Roman-numeral label and a title in small caps (e.g., "I   ALTERNATIVES",
  "II   STATES OF THE WORLD")
- Within each panel, methodological choices are arranged as small rectangular chips with thin borders,
  arranged in rows, connected by implicit contrast (left = constrained/assumed, right =
  exploratory/searched)
- No directional arrows between panels — the four axes are independent
- Caption appears below, italicized, explaining each panel in numbered list form

**Elements to reuse for Figure 1 rewrite:**
- Small-caps panel headers with Roman numerals
- Rectangular chip layout for methodological sub-choices within each axis
- Clear visual distinction between the "assumed/prespecified" and "searched/discovered" sub-choices
  within each axis
- Pale monochromatic palette to keep the taxonomy readable

**What this figure does NOT show:** It does not show workflow flow — it is a comparison matrix, not
a process diagram. Do not conflate it with a flow diagram.

---

### Figure 2: `moallemi_eta_2020_exploratory_modeling_Box1.jpg`

**Source:** Moallemi et al. (2020a), Box 1 (Global Environmental Change).

**Visual grammar:**
- Two side-by-side text-diagram blocks labeled (a) and (b)
- Each block has: problem type statement (bold), decision specification, uncertainty characterization,
  outcome implication, and example (all as labeled bullet points)
- Small conceptual diagram embedded in each block: shows three overlapping regions labeled
  "Decision space", "Uncertainty space", "Outcome space" with an icon for the pre-specified alternatives
  (labeled A, B, C) and sampling arrows from the uncertainty space
- For (a) Design of experiments: sampling arrows go into the uncertainty space with no directionality
  constraint — a cloud of dots
- For (b) Stress-testing: sampling arrows plus a subspace-partitioning visualization (shaded region
  identifying failure zone)

**Elements to reuse:**
- The three-space conceptual diagram (Decision space / Uncertainty space / Outcome space) is canonical
  vocabulary for the field; our Figure 1 should use these same three labeled regions
- The visual of a cloud of sampled points vs. structured coverage in the uncertainty/outcome space is
  directly useful for contrasting MOEA-FIND with naive sampling

---

### Figure 3: `moallemi_eta_2020_exploratory_modeling_Box2.jpg`

**Source:** Moallemi et al. (2020a), Box 2 (Global Environmental Change).

**Visual grammar:**
- Three more blocks: (c) Worst-case scenario discovery, (d) Many-objective optimization, (e)
  Many-objective robust optimization
- Same structured text format as Box 1
- Small conceptual diagrams:
  - (c): Search arrows pointing toward worst-outcome region of outcome space
  - (d): Optimization arrows in decision space with Pareto front visualized in outcome space; pre-specified
    scenarios shown as points in uncertainty space
  - (e): Combined directed search in decision space with sampling in uncertainty space
- A "Directed search" label distinguishes these from the "Open exploration" of Box 1

**Elements to reuse:**
- The "Directed search" vs. "Open exploration" framing
- For our Figure 1: MOEA-FIND is a "Directed search in outcome space" — the arrow directionality should
  point *into* the hazard outcome space, not into the parameter uncertainty space as in (c)–(e)

---

### Figure 4: `moallemi_etal_2020_robustness_framework.jpg`

**Source:** Moallemi et al. (2020b), Figure 2 (Environmental Modelling and Software). Adapted from
Kwakkel and Haasnoot (2019).

**Visual grammar:**
- Five rectangular boxes connected by gray directional arrows; feedback loops shown as return arrows
- Left side: two stacked boxes feeding into "Model simulations":
  - "Generation of decisions" (with sub-labels: Search — Many objective optimisation, Global or local
    sampling | Pre-specified — Expert opinion, standardised)
  - "Generation of scenarios" (same sub-label structure)
  - Each of these has a labeled output arrow: "Decision space" and "Uncertainty space" respectively
- Center: "Model simulations" box → labeled output arrow "Outcome space" → right
- Right: "Measurement of performance" box with three sub-labels (Regret, Satisficing, Descriptive stats)
- Bottom: "Vulnerability analysis" box with sub-labels (Subspace partitioning, Scenario discovery and
  adaptation tipping points | Sensitivity analysis: Factor fixing and factor prioritisation)
- Feedback: "Vulnerability analysis" has upward arrows going back to both "Generation of" boxes and to
  "Measurement of performance"

**This is the most important figure for our Figure 1 redesign.** The grammar:
- Left-to-right flow: generate inputs → simulate → measure performance
- Feedback from vulnerability analysis closes the loop
- Generation of scenarios is a distinct named box equal in status to generation of decisions

**What to modify for our Figure 1:**
- Within the "Generation of scenarios" box, distinguish sub-rows: (1) parameter-space sampling
  (conventional, shown as currently labeled), (2) MOEA-FIND — directed search in hazard-outcome space
  using Borg MOEA + Kirsch-Nowak
- Highlight that MOEA-FIND's search operates not in the "Uncertainty space" box but across the mapping
  from uncertainty space through Model simulations to Outcome space — it closes a feedback loop earlier
  than vulnerability analysis
- Keep the same overall box-and-arrow grammar and arrow style (thick gray)

---

## 4. Key Framing Nuance — MOEA-FIND and Parameter Space

**This section is load-bearing for the rewrite and must not be simplified away.**

MOEA-FIND does not replace parameter-space sampling. The relationship is nested and emergent:

1. Kirsch-Nowak generates synthetic streamflow traces by sampling *block bootstrap indices* from
   historical data. These bootstrap indices constitute the parameter space (or decision variable space
   in the optimization formulation). A given parameter vector (set of bootstrap indices) produces one
   synthetic trace, and that trace has a set of drought characteristics (duration, average severity,
   Manhattan norm from anti-ideal) which place it at a point in hazard-outcome space.

2. MOEA-FIND's Borg MOEA searches over the parameter space (bootstrap indices) to find combinations
   that produce drought characteristics with near-uniform coverage of the reachable hazard-outcome space.
   The Pareto front — in hazard space — is the structured ensemble.

3. **The hazard space is therefore EMERGENT from the parameter space.** MOEA-FIND does not bypass
   parameters; it inverts the usual relationship. Instead of: sample parameters → accept whatever hazard
   space coverage results, MOEA-FIND does: define target coverage in hazard space → search parameter
   space to achieve it.

4. **What MOEA-FIND does NOT do:** It does not sample the full joint distribution over generator
   parameters (e.g., the covariance structure used by Kirsch-Nowak, or choice of historical record
   length). A single MOEA-FIND run holds these generator parameters fixed and produces structured hazard
   coverage conditional on those parameters.

5. **Natural nested extension (Discussion section, not claimed as demonstrated):** A nested design
   would place MOEA-FIND inside an outer loop that samples generator parameters (e.g., LHS over
   Kirsch-Nowak parameters, or over historical conditioning windows). Each outer sample would produce
   one MOEA-FIND ensemble with structured hazard coverage; the collection of these ensembles would
   provide both structured hazard coverage and uncertainty quantification over generator parameters.
   This addresses the stochastic-convergence concern: without the outer loop, one MOEA-FIND run
   might find a systematically less hazardous region due to the particular historical conditioning
   of the generator, not due to anything about the decision alternative.

6. **Correct framing for the Introduction:** MOEA-FIND is a method for structured sampling of the
   emergent hazard space, rather than relying on very large parameter ensembles to achieve coverage by
   brute force. It is complementary to — not a replacement for — parameter-space uncertainty
   quantification. The appropriate Intro language is: "MOEA-FIND operates at the scenario-generation
   stage of robustness analysis [following Moallemi 2020b, Stage II, fork 3.2.2], generating states of
   the world with structured coverage of drought characteristic space. Unlike conventional parameter
   sampling, which achieves coverage of the hazard space only implicitly through large ensemble size,
   MOEA-FIND uses directed multi-objective search to ensure near-uniform representation of the
   reachable hazard space."

**Do NOT write:** "MOEA-FIND replaces parameter sampling" or "MOEA-FIND eliminates the need for
large parameter ensembles." These framings are wrong and will invite justified rejection.

---

## 5. Recommended Citation Anchors by Intro Paragraph

### P1 — Robustness problem: scenario definition is load-bearing

**Primary Herman anchors:**
- Abstract: "We propose a taxonomy of robustness frameworks... underscoring the importance of an
  informed definition of robustness."
- Intro: "the methodological choices in the taxonomy lead to the selection of substantially different
  planning alternatives"
- Conclusions: "methodological choices in the robustness analysis may lead to significantly different
  recommendations from the Pareto set"

**Primary Moallemi (2020b) anchors:**
- Abstract: "These methodological choices... can lead to fundamentally different outcomes for the
  systems of focus"
- Path dependency (Section 2.3): "final decision recommendations [are] highly sensitive to the
  methodological path taken, or path dependent"

**Use:** These passages support the claim that scenario definition choices are consequential, not
merely technical.

---

### P2 — Dominant practices: parameter-space sampling and storyline approaches

**Herman anchors:**
- Axis II discussion: "Decision Scaling, RDM, and MORDM typically perform an exploratory global
  sample of uncertain factors over plausible ranges"
- States of the world sub-axis: "Key Factors Assumed | Key Factors Discovered" and
  "Prespecified Scenarios | Design of experiments"

**Moallemi (2020a) anchors:**
- Box 1a (Design of experiments): "random values from the uncertainty space are systematically sampled"
- Box 1b (Stress-testing): same parameter-sampling characterization
- Section 2: "Different exploratory modeling approaches share the core idea of systematically
  analyzing the implications of decision and uncertainty spaces in the outcome space"

**Moallemi (2020b) anchors:**
- Section 3.2.2: the two main modes for framing scenarios are "pre-specified standardised
  (and deterministic) frames of the future" and "globally sampled stochastic analysis"

These passages characterize parameter-space sampling and pre-specified storylines as the two dominant
practices — exactly the P2 claim the Introduction needs to make.

---

### P3 — Hazard space as explicit target; prior work

**Gap note:** Neither Herman (2015) nor Moallemi (2020a, 2020b) discuss hazard-outcome space as a
sampling target. Their scenario generation is entirely framed around *parameter* space. The
hazard-outcome space as an explicit target for scenario generation is the specific conceptual
contribution of Borgomeo et al. (2015), Zaniolo et al. (2024), and the MOEA-FIND paper.

**Recommended approach:** Do not cite Herman or Moallemi for the hazard-space framing. The Introduction
should say something like: "While [Herman 2015] and [Moallemi 2020] identify scenario generation as a
load-bearing methodological choice, neither framework provides methods for targeting the coverage of
drought hazard characteristics directly. Borgomeo et al. (2015) and Zaniolo et al. (2024) introduced
scenario-generation methods that operate in drought characteristic space..." then position MOEA-FIND
as extending this line.

---

### P4 — MOEA-FIND positioning in Herman and Moallemi frameworks

**Herman (2015):**
- Axis II (States of the World): MOEA-FIND is a search-based method for generating states of the world,
  filling the "Design of experiments → Search" sub-choice that Herman only partially maps. No direct
  passage supports this; it must be argued by analogy to the axis structure.
- The key passage: "This sampling is typically performed over noninformative priors, reflecting the
  exploratory nature of deep uncertainty analyses" — MOEA-FIND replaces noninformative parameter
  sampling with directed hazard-space search.

**Moallemi (2020b) Section 3.2.2:**
- "Framing future scenarios (i.e., states of the world) over which decisions are to be evaluated" is
  the named stage where MOEA-FIND operates.
- The two existing modes (pre-specified and globally sampled) do not include directed hazard-space
  search — MOEA-FIND is the third mode for this fork.

**Recommended Introduction claim:** "In the taxonomy of [Herman 2015] and the decision support
framework of [Moallemi 2020b], scenario generation constitutes a distinct methodological stage —
Axis II (States of the World) and Stage II fork 3.2.2 (Framing future scenarios), respectively. Prior
methods at this stage sample the parameter uncertainty space. MOEA-FIND introduces directed search
in the emergent hazard-outcome space as an alternative mode of scenario generation within this stage."

---

### P5 — Scope and critique defusal

**Herman anchors:**
- "analyses must solidify the link to top-down projections such as climate forecasts to better
  constrain the range of plausible future states of the world, including estimates of likelihoods
  where appropriate" — Herman himself notes this open question; MOEA-FIND does not resolve it and
  should not claim to
- The paper's overall mantra: "key uncertainties should be discovered through sensitivity analysis
  rather than assumed" — MOEA-FIND discovers hazard coverage rather than assuming it

**Moallemi (2020b) anchors:**
- The path-dependency framing (Section 2.3): MOEA-FIND is a new path at the scenario-generation fork,
  not a claim that all other paths are wrong
- Section 4 (challenges): "computational limitations in working with large-scale assessment models" —
  MOEA-FIND has the same computational challenge; the HPC context should be noted

---

## 6. Vocabulary Checklist

The following phrases appear verbatim in Herman (2015) and/or Moallemi (2020) and should be used
exactly in the Introduction to signal fluency with this literature:

**From Herman (2015):**
1. "states of the world" — use throughout (not "scenarios" alone, not "futures")
2. "the selection and sampling of the uncertain factors which make up plausible future states of the world" — this is the exact framing of Axis II; quote or closely paraphrase
3. "sampling of states of the world" — the axis name Herman uses
4. "factor mapping" — sensitivity analysis that identifies which uncertain factors cause failure
5. "factor prioritization" — sensitivity analysis that ranks all uncertain factors
6. "domain criterion" — the satisficing robustness measure (fraction of states where requirements are met)
7. "a posteriori decision support" — discovering alternatives/scenarios before imposing preferences
8. "generate-first-choose-later" — the GFCL mantra
9. "deep uncertainty" / "deeply uncertain" — use consistently, not "high uncertainty" or "large uncertainty"
10. "scenario-neutral" — Prudhomme et al. (2010), cited in Herman — relevant if discussing bottom-up framing

**From Moallemi (2020a, Global Environmental Change):**
11. "uncertainty characterization" — the component describing how scenarios are generated
12. "decision specification" — the component describing how alternatives are generated
13. "outcome space" — where model simulation results live; distinguish from uncertainty space and decision space
14. "open exploration" vs. "directed search" — these are Moallemi's exact terms for the two modes

**From Moallemi (2020b, Environmental Modelling and Software):**
15. "framing future scenarios" — the exact name of Stage II fork 3.2.2 where MOEA-FIND operates
16. "generation of scenarios (i.e., states of the world)" — the exact axis name in the nested taxonomy
17. "path dependency" — the sensitivity of final inferences to methodological choices at each fork
18. "decision forks" — the named methodological choice points in the decision support process

---

*Document prepared 2026-04-15. No Introduction text drafted. Existing manuscript files not modified.*
