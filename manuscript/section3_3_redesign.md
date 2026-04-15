# Section 3.3 Redesign: Scenario Discovery with Operational Failure Labels

*Written 2026-04-14. Purpose: replace the circular severity/duration threshold label
identified in HC-2 (reviewer_critiques.md) with an operationally grounded failure label
derived from Pywr-DRB simulation. No runs are performed here — this is a paper-design
document specifying inputs, configuration, metrics, figure layout, and the claim the
redesigned section supports.*

---

## 1. The problem with the current design

Section 3.3 currently defines a binary failure label as the joint condition that
mean drought severity $D_1$ and mean drought duration $D_2$ both exceed their
respective historical 80th percentiles, and then trains a gradient boosted tree
classifier on $D_1$ and $D_2$ as features against that label. The circularity is
exact: the failure threshold is defined on the same coordinate system that MOEA-FIND
optimizes over, so any space-filling sample in $(D_1, D_2)$ will produce a
better-calibrated classifier of a threshold in $(D_1, D_2)$ than a density-weighted
library subsample will. This is a property of the sampling design, not a property of
MOEA-FIND specifically.

The redesign replaces the threshold label with a binary outcome derived from
Pywr-DRB simulation: a drought scenario causes a system failure if the New York City
water supply system, driven with that scenario's Cannonsville inflows and a reference
set of other-site inflows, falls below an operationally defined reliability or
vulnerability threshold. The failure label is computed by a separate model from the
drought characteristics the optimizer targeted, so its relationship to the
$(D_1, D_2, D_3)$ coordinates is an empirical finding rather than a tautology.

---

## 2. Experiment inputs

### 2.1 Synthetic trace ensembles

Two ensembles of equal size $N$ enter Pywr-DRB simulation. $N$ is determined by the
MOEA-FIND archive size at convergence (anticipated 150–400 members for a
three-characteristic, single-site run at the target epsilon values). Both ensembles
share the same $N$.

**Ensemble A — MOEA-FIND archive.** The $N$ Pareto-optimal synthetic Cannonsville
inflow traces returned by the convergence of the MOEA-FIND optimization run
(Section 3.2). Each member is a monthly streamflow time series whose drought
characteristics $(D_1, D_2, D_3)$ cover the feasible drought hazard region with
near-uniform spacing under the epsilon-dominance archive.

**Ensemble B — Library LHS subsample.** $N$ traces drawn from the 10,000-trace
Kirsch-Nowak library (Section 3.2) using conditioned Latin Hypercube Sampling
(cLHS, Minasny and McBratney, 2006) applied in the generator decision-variable
space, following the approach of Bonham et al. (2024). Subsampling in input
(decision-variable) space rather than hazard space is the Bonham et al. design;
this is the appropriate baseline because it represents the best-practice alternative
for input-space ensemble reduction. The $N$ members of Ensemble B cover the
input space with near-uniform spacing but produce a non-uniform projection onto
$(D_1, D_2, D_3)$ hazard space, as the generator map is nonlinear.

### 2.2 Pywr-DRB configuration

**Inflow assembly.** Each ensemble member provides Cannonsville inflows for the
simulation period (matched to the trace length of the generated records). Inflows at
all other Pywr-DRB sites (Pepacton, Neversink, Mongaup, and main-stem lateral
inflows) are drawn from historical resampling using the same Kirsch-Nowak generator
with no drought steering, so that only the Cannonsville inflow varies across
simulation runs. This isolates the effect of Cannonsville drought characteristics on
system performance and avoids confounding from correlated multi-site drought events
at this stage. Multi-site correlation is addressed in the companion multi-site
application paper.

**Initial conditions.** All runs begin from the same initial reservoir storage
state, set to the median Cannonsville storage observed on 1 April across the
historical record (the beginning of the primary drought-season build-up period for
the NYC Delaware system). Using a fixed starting state removes initial-condition
variance as a confounding factor in the failure label.

**Operating rules.** Standard Flexible Flow Management Program (FFMP) rules as
implemented in the current Pywr-DRB model (Amestoy et al., in prep) are used
without modification. No policy optimization is performed within Section 3.3; the
experiment evaluates the performance of the existing operating policy under each
generated drought scenario.

**Demand.** A flat demand scenario equal to the average historical NYC Delaware
Basin per capita demand is applied, following the standard Pywr-DRB single-scenario
convention for baseline diagnostic runs.

**Simulation period.** Matched to the trace length of the generated records. For
the 30-year trace length used in the Cannonsville demonstration, each Pywr-DRB
run covers 30 years of monthly time steps (360 time steps per run).

---

## 3. Operational failure metrics

Three Hashimoto et al. (1982) metrics are computed from the monthly storage and
shortage time series output by each Pywr-DRB simulation run:

**Supply reliability** $\lambda$: the fraction of months in which the NYC combined
storage across the Cannonsville, Pepacton, and Neversink reservoirs remains above the
Level 1 FFMP drought trigger threshold. A month is a supply failure if combined
storage falls below that threshold. $\lambda = 1 - (\text{failure months}) / 360$.

**Vulnerability** $\nu$: the maximum single-month shortfall in deliverable supply
relative to demand, normalized by monthly demand. $\nu = \max_t \bigl[
(\text{demand}_t - \text{deliverable}_t) / \text{demand}_t \bigr]$, where deliverable
supply is constrained by FFMP-governed release rules. Months with no shortage
contribute zero to the maximum.

**Flow target reliability** $\mu$: the fraction of months in which the Montague
and Trenton non-NYC flow targets specified by the Delaware River Compact are both
met. A month is a flow-target failure if either target is missed after FFMP
diversions. $\mu = 1 - (\text{compact-failure months}) / 360$.

**Binary failure label** $y \in \{0, 1\}$: a simulation run is classified as a
system failure ($y = 1$) if supply reliability falls below the threshold
$\lambda < \lambda^*$ or vulnerability exceeds the threshold $\nu > \nu^*$. The
thresholds $\lambda^*$ and $\nu^*$ are calibrated against the historical Cannonsville
record: $\lambda^* = 0.95$ (the historical observed reliability over the 1945–2023
period) and $\nu^* = 0.15$ (the 90th percentile of monthly shortfall severity in
the historical record). Sensitivity to these thresholds is reported in Supporting
Information. The flow target reliability $\mu$ enters as a secondary label in a
two-label sensitivity analysis but is not used for the primary binary classification.

---

## 4. Classifier and comparison

**Classifier.** A gradient boosted tree classifier (Friedman, 2001; Chen and
Guestrin, 2016) is trained on each ensemble separately to predict the binary
failure label $y$ from the three drought characteristic features
$(D_1, D_2, D_3)$. Hyperparameters (maximum depth, learning rate, number of
estimators) are tuned by five-fold cross-validation within the training set and
are held fixed across both ensembles to ensure a controlled comparison.

**What is being compared.** For both ensemble A (MOEA-FIND) and ensemble B
(library LHS subsample), the classifier is trained on the same $N$ feature-label
pairs $\{(D_1^{(i)}, D_2^{(i)}, D_3^{(i)}), y^{(i)}\}_{i=1}^{N}$. The test set
is a held-out random subsample of 2,000 traces from the full 10,000-trace library,
with operational labels obtained from an independent set of Pywr-DRB simulation runs.
The test set is identical for both classifiers.

**Comparison metrics.**
- Area under the ROC curve (AUC) on the test set.
- Brier score (mean squared probability calibration error) on the test set.
- Decision boundary visualization: the classifier's predicted failure probability
  as a function of $(D_1, D_2)$ at the median value of $D_3$, displayed as a
  filled contour overlaid on the ensemble scatter plots.
- Near-boundary coverage: the fraction of training points within a narrow band
  of the true failure boundary (estimated from the test set classifier, which
  has $20\times$ more training data). A higher near-boundary fraction indicates
  that the training set has more points informative for locating the failure
  threshold.

---

## 5. Figure layout

**Figure 3.3** (four-panel, landscape orientation, approximately 160 mm wide).

**Panel (a).** Scatter plot of Ensemble A (MOEA-FIND, $N$ points) in the
$(D_1, D_2)$ plane, with $D_3$ encoded by symbol size. Points are colored by
operational failure label: red for $y=1$, blue for $y=0$. The convex hull of
the feasible drought hazard region estimated from the 10,000-trace library is
shown as a light gray polygon. The MOEA-FIND archive is expected to fill this
region with near-uniform spacing, with both failure and non-failure points
distributed across the region including the interior.

**Panel (b).** Same axes and color conventions as panel (a), showing Ensemble B
(library LHS subsample, $N$ points). The library LHS subsample is expected to show
density weighting toward moderate drought conditions, with sparse coverage near the
tails of the feasible region where the failure boundary is likely to lie.

**Panel (c).** Decision boundary comparison. The predicted failure probability
from the MOEA-FIND-trained classifier and from the LHS-trained classifier are both
shown as filled contours in the $(D_1, D_2)$ plane at median $D_3$. The true
failure boundary estimated from the large test set is overlaid as a solid contour
line. A well-calibrated classifier closely follows the true boundary; a
poorly-calibrated classifier shows a diffuse or displaced boundary. The expected
finding is that the MOEA-FIND classifier tracks the true boundary more closely
because its training data cover both sides of the boundary at multiple points.

**Panel (d).** Bar chart or strip chart comparing four scalar metrics across the
two classifiers (MOEA-FIND vs. LHS subsample): AUC, Brier score, near-boundary
training fraction, and failure-region boundary length recovered. Error bars from
five-fold cross-validation.

---

## 6. Claim the redesigned section supports

The primary claim is: synthetic trace ensembles with structured coverage of the
drought hazard space produce more accurate identification of the operational failure
boundary in that space than equal-size ensembles sampled for coverage of the
generator input space, because structured hazard-space coverage ensures training
data on both sides of the failure boundary across the full range of drought
characteristics.

This claim is empirically testable and non-trivial: if the operational failure
boundary in $(D_1, D_2, D_3)$ space is simple (e.g., a halfspace at high severity
and duration), an LHS subsample from the tail of the input distribution may perform
comparably. If the boundary is curved or has interior regions of non-failure
surrounded by failure (e.g., because moderate droughts during certain seasonal
windows are more damaging than severe off-season droughts), MOEA-FIND's interior
coverage will be essential for boundary recovery. The finding is informative
regardless of direction.

The redesigned section explicitly does not claim that MOEA-FIND discovers more
operational failures overall than library subsampling. It claims that MOEA-FIND
produces a better-calibrated map of which drought conditions cause failure. This
distinction is important: finding more failure scenarios is not the goal; locating
the boundary of the failure region in drought hazard space is.

---

## 7. Relationship to HC-2 and the complementary framing

The revised Section 3.3 addresses HC-2 (circularity) by ensuring the failure label
is derived from a model (Pywr-DRB) whose outputs are not defined in terms of the
drought characteristic coordinates. The Pywr-DRB model maps inflow time series to
reservoir storage and shortage time series through physical routing and operating
rules; the relationship between $(D_1, D_2, D_3)$ and the Hashimoto reliability and
vulnerability metrics is a property of the basin hydrology, reservoir geometry,
operating rules, and demand structure, not of the drought metric definitions.

The experiment also reinforces the complementary framing established in the
Introduction. Section 3.3 uses Pywr-DRB with a fixed operating policy, not a
policy search. The role of the MOEA-FIND ensemble in this context is to supply
the scenario basis for vulnerability discovery, not to replace the simulation model
or the policy optimization. The workflow is: (MOEA-FIND generates structured
ensemble) + (Pywr-DRB evaluates each member) + (scenario discovery in hazard
space identifies the failure boundary). MOEA-FIND contributes the scenario design
step; the other two steps are unchanged from standard practice.

---

## 8. Pending dependencies and sequencing

This experiment depends on two prior completed components:

1. **MOEA-FIND archive (Section 3.2, Phase beta HPC):** The $N$ archive members
   must be available as monthly inflow time series files compatible with Pywr-DRB
   input format.

2. **Pywr-DRB simulation runs:** Each of the $N_A + N_B + 2000$ traces requires one
   Pywr-DRB simulation run. At approximately 150 seconds per 30-year run, the total
   compute cost is approximately $(N_A + N_B + 2000) \times 150$ seconds. For
   $N = 300$ and the 2,000-trace test set, total wall-clock time is approximately
   390,000 seconds (108 hours) in serial, or under 4 hours with 32-core parallelism.
   This is within the HPC Phase gamma budget.

3. **Library operational labels (test set):** The 2,000-trace test set requires
   independent Pywr-DRB runs from the 10,000-trace library. These can be batched
   with Phase beta HPC runs.
