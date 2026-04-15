# HPC Deployment Status — MOEA-FIND Sections 3.2 and 3.3

*Created: 2026-04-14 by autonomous HPC deployment session (eager-dirac worktree).*
*Purpose: survey what exists vs. what was created in this session, and flag TODOs for Trev.*

---

## Summary

The pipeline for Section 3.2 (Cannonsville coverage comparison) and Section 3.3
(Pywr-DRB reliability/vulnerability labeling) is now **deployable** in the sense
that all scripts are written and the preflight test is ready to run. The scripts
depend on design decisions listed under TODOs below that must be resolved before
submitting production HPC jobs.

**Run the preflight test first:**
```bash
python scripts/hpc/preflight_test.py --skip-pywrdrb   # without Pywr-DRB
python scripts/hpc/preflight_test.py                   # with Pywr-DRB installed
```

If preflight passes, the submission sequence is:
```bash
# Section 3.2 — can run concurrently
sbatch scripts/hpc/section32_cannonsville_moea.slurm
sbatch scripts/hpc/section32_cannonsville_library.slurm

# Section 3.3 — requires both section32 jobs to complete first
sbatch scripts/hpc/section33_pywr_drb_batch.slurm
```

---

## What Existed Before This Session

| Item | Location | Status |
|------|----------|--------|
| KirschBorgWrapper | `src/kirsch_wrapper.py` | Complete |
| LibraryGenerator | `src/library.py` | Complete |
| SSI objectives + Manhattan norm | `src/objectives.py` | Complete |
| Coverage metrics (NN-CV, L2*) | `src/analysis.py` | Complete |
| Shared experiment utilities | `src/experiment_utils.py` | Complete |
| USGS data loader (Cannonsville) | `src/data.py` | Complete |
| Kirsch single-site MOEA driver | `scripts/04_kirsch_single_site.py` | Complete (but not config-driven, saves objs only) |
| Library generation (MPI) | `scripts/05_kirsch_library_build.py` | Complete (saves characteristics, NOT raw flows) |
| Library LHS subsample | `scripts/06_library_subsample_baseline.py` | Complete (hazard-space LHS only) |
| Pywr-DRB policy re-eval | `scripts/09_drb_policy_reeval.py` | Scaffolded stub only |
| SLURM common helpers | `scripts/_common.sh` | Complete |
| Existing SLURM scripts (01–10) | `scripts/NN_*.slurm` | Complete for analytic/Kirsch experiments |

---

## What Was Created in This Session

| Item | Location | Notes |
|------|----------|-------|
| Section 3.2 config | `configs/section32_cannonsville.yaml` | All tunables in one place |
| Section 3.3 config | `configs/section33_pywr_drb.yaml` | Pywr-DRB settings, thresholds |
| Section 3.2 MOEA driver | `scripts/hpc/section32_cannonsville_moea.py` | Config-driven, saves flow arrays for Pywr-DRB |
| Section 3.2 MOEA SLURM | `scripts/hpc/section32_cannonsville_moea.slurm` | EpsNSGAII serial; MM Borg TODO |
| Section 3.2 library SLURM | `scripts/hpc/section32_cannonsville_library.slurm` | Calls scripts 05 + 06; saves npz |
| Section 3.3 Pywr-DRB batch | `scripts/hpc/section33_pywr_drb_batch.py` | mpi4py workers, inflow injection, metric extraction |
| Section 3.3 Pywr-DRB SLURM | `scripts/hpc/section33_pywr_drb_batch.slurm` | 32 cores, 8 h allocation |
| Pre-flight test | `scripts/hpc/preflight_test.py` | 10-step pipeline check, login-node safe |
| HPC conda environment | `environment_hpc.yml` | Pins for Python 3.11, MPI, Pywr-DRB |

---

## TODOs for Trev (in priority order)

These are either methodology decisions or system-specific configuration items.
None of these were changed by this session.

### T-01 (BLOCKING for Section 3.3): Pywr-DRB inflow injection approach
**File:** `scripts/hpc/section33_pywr_drb_batch.py` — `build_custom_inflow_dir()` and `run_pywrdrb_simulation()`

The batch wrapper creates a temporary flow directory per trace and registers it
with pywrdrb's PathNavigator via `pn.sc.add(f"flows/{inflow_type_name}", ...)`.
This pattern is used in `model_builder.py` for the ML plugin. Confirm it works
for the base model builder (i.e., that all required files are found in the temp dir).

Required files in the temp directory:
- `catchment_inflow_mgd.csv` (Cannonsville column replaced ✓)
- `gage_flow_mgd.csv` (copied unchanged — **TREV: confirm Cannonsville gage_flow should also be replaced**)
- `predicted_inflows_mgd.csv` (copied unchanged — affects STARFIT)
- Possibly: `diversion_nyc_extrapolated_mgd.csv`, `diversion_nj_extrapolated_mgd.csv`

**Action:** Run preflight step 9 on HPC, then test with 1 live simulation
(`python scripts/hpc/section33_pywr_drb_batch.py --mode moea --max-traces 1`).

### T-02 (BLOCKING for Section 3.3): Monthly→daily disaggregation
**File:** `configs/section33_pywr_drb.yaml` → `pywr_drb.monthly_to_daily`

Current placeholder: `"uniform"` — each day in month receives the monthly mean.
**section3_3_redesign.md** does not specify a disaggregation method. Options:
- `"uniform"`: Simplest. Removes daily variability. Acceptable for monthly Pywr-DRB analysis.
- `"ratio"`: Multiply historical daily pattern by the ratio of synthetic to historical monthly mean. Preserves day-of-month variability but mixes historical and synthetic signals.

**Action:** Choose method before final analysis. Change `monthly_to_daily` in the config.

### T-03 (BLOCKING for Section 3.3): Confirm Pywr-DRB simulation period
**File:** `configs/section33_pywr_drb.yaml` → `pywr_drb.sim_start_date`, `sim_end_date`

Current: `1980-10-01` to `2010-09-30` (30 water years within pub_nhmv10_BC_withObsScaled).
The base inflow type covers 1945–2022. Any 30-year window within that range works.
The choice affects what "other-site" flows the non-Cannonsville catchments receive.

**Action:** Confirm or adjust the 30-year window.

### T-04 (BLOCKING for Section 3.3): Compute initial_volume_frac from historical record
**File:** `configs/section33_pywr_drb.yaml` → `pywr_drb.initial_volume_frac`

Current placeholder: `0.8`.
**section3_3_redesign.md §2.2:** "median Cannonsville storage observed on 1 April."

To compute:
```python
import pywrdrb
# Run a full historical simulation, load res_storage, extract
# Cannonsville values on April 1 of each year, compute median / max_storage.
```

**Action:** Run a historical Pywr-DRB simulation, extract April 1 storage fraction.

### T-05 (BLOCKING for Section 3.3): Confirm Pywr-DRB output column mapping
**File:** `configs/section33_pywr_drb.yaml` → `metrics.*`
**File:** `scripts/hpc/section33_pywr_drb_batch.py` — `extract_metrics()`

The metric extraction uses these Pywr-DRB output keys:
- `res_storage` with columns `cannonsville`, `pepacton`, `neversink`
- `res_level` with column `res_level` (FFMP drought level — **confirm encoding: 0=normal, 1=L1, 2=L2, ...?**)
- `major_flow` with columns `delMontague`, `delTrenton`
- `mrf_target` with columns matching `major_flow`
- `nyc_deliveries` for ν computation (**TREV: confirm this key exists in current pywrdrb schema**)

**Action:** Run the historical simulation used in T-04 with the relevant results_sets
and print `data.res_level[label][0].columns` etc. to confirm mappings.

### T-06 (METHODOLOGY — Trev): Third drought objective for k=3 Pareto front
**File:** `configs/section32_cannonsville.yaml` → `moea.objective_keys`

**section3_3_redesign.md** refers to three characteristics (D1, D2, D3) but does
not explicitly define D3. Current MOEA driver uses k=2 (`mean_duration`, `mean_avg_severity`).
To run k=3, uncomment the third entry in `objective_keys` in the config (e.g., `frequency`).
**No code change needed — config change only.**

**Action:** Choose D3 and update `objective_keys` in `configs/section32_cannonsville.yaml`.

### T-07 (METHODOLOGY — Trev): Baseline subsampling method (Bonham et al. DV-space)
**File:** `configs/section32_cannonsville.yaml` → `library.subsample_method`
**File:** `scripts/hpc/section32_cannonsville_library.slurm`

Current: hazard-space LHS (nearest-neighbor matching in drought characteristic space).
**section3_3_redesign.md §2.1:** baseline should use cLHS in decision-variable space
(Bonham et al., 2024). This changes the subsampling from `scripts/06_library_subsample_baseline.py`
to a new cLHS-in-DV-space implementation.

**Action:** Decide which baseline to use for submission and implement cLHS if chosen.

### T-08 (METHODOLOGY — Trev): Library trace raw flow saving (script 05)
**File:** `scripts/05_kirsch_library_build.py`

Script 05 saves `characteristics.json` but not the raw flow arrays.
Section 3.3 batch wrapper needs `library_traces.npz` (raw flows) to pass
the test set to Pywr-DRB. The `section32_cannonsville_library.slurm` notes
this gap with a placeholder `flows_cfs=np.zeros(...)`.

**Action:** Modify `scripts/05_kirsch_library_build.py` to also save `library_traces.npz`.
(This is an execution detail, not a methodology decision — safe to do.)

### T-09 (SYSTEM): Edit `scripts/_common.sh` for HPC cluster
**File:** `scripts/_common.sh`

Cluster-specific variables at the top:
```bash
CLUSTER_ACCOUNT="CHANGE_ME_ACCOUNT"
CLUSTER_PARTITION="CHANGE_ME_PARTITION"
CLUSTER_PYTHON_MODULE="python/3.11"
CLUSTER_MPI_MODULE="openmpi/4.1.5"
CLUSTER_VENV="${HOME}/.venvs/moea-find"
```

**Action:** Fill in account name, partition name, and module names for your HPC.

### T-10 (SYSTEM): MM Borg Python wrapper installation
**File:** `scripts/hpc/section32_cannonsville_moea.slurm`

Current MOEA driver uses EpsNSGAII (serial). MM Borg requires `borg.py` and
compiled MPI binaries on PYTHONPATH. The slurm comment notes this.

**Action:** See Feb 2025 and Aug 2025 Water Programming blog posts (noted in CLAUDE.md)
for installation instructions. Once installed, change the slurm launch to use
`moea_mpi_launch` and update the Python driver to import Borg instead of platypus.

---

## Compute Budget

| Job | Cores | Wall time | Notes |
|-----|-------|-----------|-------|
| section32 MOEA | 1–64 | 8 h | Serial EpsNSGAII; extra nodes idle unless Borg used |
| section32 library | 32 | 6 h | MPI: 10K traces / 32 workers ≈ 313 traces/worker |
| section33 Pywr-DRB | 32 | 8 h | ~2600 traces × 150 s / 32 ≈ 3.4 h; 8 h = 2× buffer |

Budget source: section3_3_redesign.md §8 (Phase gamma HPC budget).

---

## File Tree (new files from this session)

```
configs/
  section32_cannonsville.yaml
  section33_pywr_drb.yaml
scripts/hpc/
  section32_cannonsville_moea.py
  section32_cannonsville_moea.slurm
  section32_cannonsville_library.slurm
  section33_pywr_drb_batch.py
  section33_pywr_drb_batch.slurm
  preflight_test.py
manuscript/
  hpc_deployment_status.md   (this file)
environment_hpc.yml
```
