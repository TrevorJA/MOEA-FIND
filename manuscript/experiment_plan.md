# Experiment Plan — Sections 3.2 and 3.3

*Merged 2026-04-15 from `section3_3_redesign.md` and `hpc_deployment_status.md`.*

---

## Part A: Section 3.3 Redesign — Scenario Discovery with Operational Failure Labels

*Written 2026-04-14. Replaces the circular severity/duration threshold label identified
in HC-2 (reviewer_defenses.md) with an operationally grounded failure label derived from
Pywr-DRB simulation.*

### A.1 The problem with the current design

Section 3.3 currently defines a binary failure label as the joint condition that mean
drought severity $D_1$ and mean drought duration $D_2$ both exceed their respective
historical 80th percentiles, and then trains a gradient boosted tree classifier on
$D_1$ and $D_2$ as features against that label. The circularity is exact: the failure
threshold is defined on the same coordinate system that MOEA-FIND optimizes over, so
any space-filling sample in $(D_1, D_2)$ will produce a better-calibrated classifier
of a threshold in $(D_1, D_2)$ than a density-weighted library subsample will. This is
a property of the sampling design, not of MOEA-FIND specifically.

The redesign replaces the threshold label with a binary outcome derived from Pywr-DRB
simulation: a drought scenario causes a system failure if the NYC water supply system,
driven with that scenario's Cannonsville inflows and a reference set of other-site
inflows, falls below an operationally defined reliability or vulnerability threshold.
The failure label is computed by a separate model from the drought characteristics the
optimizer targeted, so its relationship to the $(D_1, D_2, D_3)$ coordinates is an
empirical finding rather than a tautology.

### A.2 Experiment inputs

**Synthetic trace ensembles.** Two ensembles of equal size $N$ enter Pywr-DRB simulation.
$N$ is determined by the MOEA-FIND archive size at convergence (anticipated 150–400 members
for a three-characteristic, single-site run at the target epsilon values).

**Ensemble A — MOEA-FIND archive.** The $N$ Pareto-optimal synthetic Cannonsville inflow
traces returned by the convergence of the MOEA-FIND optimization run (Section 3.2). Each
member is a monthly streamflow time series whose drought characteristics $(D_1, D_2, D_3)$
cover the feasible drought hazard region with near-uniform spacing under the
epsilon-dominance archive.

**Ensemble B — Library LHS subsample.** $N$ traces drawn from the 10,000-trace
Kirsch-Nowak library using conditioned Latin Hypercube Sampling (cLHS, Minasny and
McBratney, 2006) applied in the generator decision-variable space, following Bonham
et al. (2024). This is the appropriate baseline: it represents best-practice
input-space ensemble reduction. The $N$ members cover the input space with near-uniform
spacing but produce a non-uniform projection onto $(D_1, D_2, D_3)$ hazard space.

**Pywr-DRB configuration.**
- Each ensemble member provides Cannonsville inflows; all other sites draw from historical
  resampling using the same Kirsch-Nowak generator without drought steering.
- Initial conditions: median Cannonsville storage observed on 1 April across the historical
  record (computed from a historical Pywr-DRB simulation — see T-04 below).
- Operating rules: standard FFMP rules as implemented in the current Pywr-DRB model.
- Demand: flat scenario equal to average historical NYC Delaware Basin per capita demand.
- Simulation period: matched to trace length (30 years = 360 monthly time steps per run).

### A.3 Operational failure metrics

Three Hashimoto et al. (1982) metrics:

**Supply reliability** $\lambda$: fraction of months in which NYC combined storage
(Cannonsville + Pepacton + Neversink) remains above the Level 1 FFMP drought trigger.
$\lambda = 1 - (\text{failure months}) / 360$.

**Vulnerability** $\nu$: maximum single-month shortfall in deliverable supply relative to
demand, normalized by monthly demand. $\nu = \max_t [(\text{demand}_t - \text{deliverable}_t)
/ \text{demand}_t]$.

**Flow target reliability** $\mu$: fraction of months in which Montague and Trenton
non-NYC flow targets are both met. Secondary label only.

**Binary failure label** $y \in \{0, 1\}$: $y = 1$ if $\lambda < \lambda^*$ or
$\nu > \nu^*$. Thresholds calibrated from the historical record: $\lambda^* = 0.95$
(historical observed reliability, 1945–2023) and $\nu^* = 0.15$ (90th percentile of
monthly shortfall severity). Threshold sensitivity reported in Supporting Information.

### A.4 Classifier and comparison

A gradient boosted tree classifier (Friedman, 2001; Chen and Guestrin, 2016) is trained
on each ensemble separately to predict the binary failure label $y$ from the three drought
characteristic features $(D_1, D_2, D_3)$. Hyperparameters are tuned by five-fold
cross-validation and held fixed across both ensembles. Test set: 2,000 held-out traces from
the full 10,000-trace library with independent Pywr-DRB operational labels (identical test
set for both classifiers).

**Comparison metrics:** AUC on test set; Brier score; decision boundary visualization;
near-boundary training fraction.

### A.5 Claim the redesigned section supports

Synthetic trace ensembles with structured coverage of the drought hazard space produce
more accurate identification of the operational failure boundary than equal-size ensembles
sampled for coverage of the generator input space, because structured hazard-space coverage
ensures training data on both sides of the failure boundary across the full range of drought
characteristics. The finding is informative regardless of direction.

### A.6 Compute budget

At approximately 150 seconds per 30-year run:
$(N_A + N_B + 2000) \times 150$ s total. For $N = 300$ and the 2,000-trace test set:
approximately 108 hours serial, or under 4 hours with 32-core parallelism.

---

## Part B: HPC Deployment Status

*Created 2026-04-14 by autonomous HPC deployment session. Records the pipeline for
Sections 3.2 and 3.3 and flags TODOs for Trev.*

### B.1 Summary

The pipeline is **deployable**: all scripts are written and the preflight test is ready.
Design decisions listed under TODOs must be resolved before submitting production HPC jobs.

**Run the preflight test first:**
```bash
python scripts/hpc/preflight_test.py --skip-pywrdrb   # without Pywr-DRB
python scripts/hpc/preflight_test.py                   # with Pywr-DRB installed
```

**Submission sequence (after preflight passes):**
```bash
# Section 3.2 — can run concurrently
sbatch scripts/hpc/section32_cannonsville_moea.slurm
sbatch scripts/hpc/section32_cannonsville_library.slurm

# Section 3.3 — requires both section32 jobs to complete first
sbatch scripts/hpc/section33_pywr_drb_batch.slurm
```

### B.2 What was created in the HPC session

| Item | Location |
|---|---|
| Section 3.2 config | `configs/section32_cannonsville.yaml` |
| Section 3.3 config | `configs/section33_pywr_drb.yaml` |
| Section 3.2 MOEA driver | `scripts/hpc/section32_cannonsville_moea.py` |
| Section 3.2 MOEA SLURM | `scripts/hpc/section32_cannonsville_moea.slurm` |
| Section 3.2 library SLURM | `scripts/hpc/section32_cannonsville_library.slurm` |
| Section 3.3 Pywr-DRB batch | `scripts/hpc/section33_pywr_drb_batch.py` |
| Section 3.3 Pywr-DRB SLURM | `scripts/hpc/section33_pywr_drb_batch.slurm` |
| Pre-flight test | `scripts/hpc/preflight_test.py` |
| HPC conda environment | `environment_hpc.yml` |

### B.3 TODOs for Trev (in priority order)

**T-01 (BLOCKING for §3.3): Pywr-DRB inflow injection approach.**
`scripts/hpc/section33_pywr_drb_batch.py` — `build_custom_inflow_dir()`. Confirm the
temporary flow directory pattern (replacing `catchment_inflow_mgd.csv` Cannonsville
column) works for the base model builder. Also confirm whether `gage_flow_mgd.csv`
should be replaced for Cannonsville. Run preflight step 9, then test with 1 live simulation.

**T-02 (BLOCKING for §3.3): Monthly→daily disaggregation.**
`configs/section33_pywr_drb.yaml` → `pywr_drb.monthly_to_daily`. Current placeholder:
`"uniform"` (each day receives the monthly mean). Alternative: `"ratio"` (multiply
historical daily pattern by ratio of synthetic to historical monthly mean). Choose before
final analysis.

**T-03 (BLOCKING for §3.3): Confirm Pywr-DRB simulation period.**
Current: `1980-10-01` to `2010-09-30` (30 water years within pub_nhmv10_BC_withObsScaled).
Confirm or adjust the 30-year window.

**T-04 (BLOCKING for §3.3): Compute `initial_volume_frac` from historical record.**
Current placeholder: `0.8`. Target: median Cannonsville storage observed on 1 April
(per §A.2 above). Run a full historical Pywr-DRB simulation, extract Cannonsville
April 1 storage fraction.

**T-05 (BLOCKING for §3.3): Confirm Pywr-DRB output column mapping.**
Confirm encoding of `res_level` (0=normal, 1=L1, 2=L2, ...?), existence of
`nyc_deliveries` key, and all column names in `res_storage`, `major_flow`, `mrf_target`.
Print columns from the historical simulation used in T-04.

**T-06 (METHODOLOGY): Third drought objective for k=3 Pareto front.**
`configs/section32_cannonsville.yaml` → `moea.objective_keys`. Current: k=2
(`mean_duration`, `mean_avg_severity`). To run k=3, uncomment the third entry (e.g.,
`frequency`). No code change needed — config change only.

**T-07 (METHODOLOGY): Baseline subsampling method.**
`configs/section32_cannonsville.yaml` → `library.subsample_method`. Current: hazard-space
LHS. §A.2 above specifies cLHS in decision-variable space (Bonham et al., 2024). Decide
which baseline to use for submission and implement cLHS in DV-space if chosen.

**T-08 (EXECUTION): Library trace raw flow saving.**
`scripts/05_kirsch_library_build.py` saves `characteristics.json` but not raw flow
arrays. §3.3 batch wrapper needs `library_traces.npz`. Modify script 05 to also save
`library_traces.npz`. Safe to do immediately.

**T-09 (SYSTEM): Edit `scripts/_common.sh` for HPC cluster.**
Fill in `CLUSTER_ACCOUNT`, `CLUSTER_PARTITION`, `CLUSTER_PYTHON_MODULE`,
`CLUSTER_MPI_MODULE`, `CLUSTER_VENV` for your HPC.

**T-10 (SYSTEM): MM Borg Python wrapper installation.**
Current MOEA driver uses EpsNSGAII (serial). Install MM Borg per Feb 2025 and Aug 2025
Water Programming posts (noted in CLAUDE.md). Change slurm launch to use
`moea_mpi_launch` and update Python driver to import Borg instead of platypus.

### B.4 Compute budget

| Job | Cores | Wall time | Notes |
|---|---|---|---|
| section32 MOEA | 1–64 | 8 h | Serial EpsNSGAII; extra nodes idle unless Borg used |
| section32 library | 32 | 6 h | MPI: 10K traces / 32 workers ≈ 313 traces/worker |
| section33 Pywr-DRB | 32 | 8 h | ~2600 traces × 150 s / 32 ≈ 3.4 h; 8 h = 2× buffer |
