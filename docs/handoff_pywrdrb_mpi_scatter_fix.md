# Handoff: Fix MPI scatter overflow in `PredictedInflowEnsemblePreprocessor.load()`

**Audience:** a fresh Claude Opus agent working in the Pywr-DRB repo (`/home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/Pywr-DRB/`).
**Goal:** replace the rank-0-scatter pattern in `src/pywrdrb/pre/predict_inflows.py` with a per-rank HDF5 read so MPI runs scale past ~16 ranks without hitting pickle-scatter buffer limits.
**Scope:** one method (`PredictedInflowEnsemblePreprocessor.load()`), ~25 lines of code, plus a smoke test. No public API change.

---

## 1. Problem

Running MOEA-FIND's step 09 (`workflows/slurm/09_drb_policy_reeval.slurm` — Stage 2 `prep`) with 120 MPI ranks over 1446 ensemble realizations fails with:

```
  File ".../pywrdrb/pre/predict_inflows.py", line 531, in load
    self.realization_data = self.comm.scatter(slices, root=0)
mpi4py.MPI.Exception: MPI_ERR_ARG: invalid argument of some other kind
```

- **Failing job:** `moea_09_reeval` 209166, Hopper, OpenMPI over InfiniBand, 3 nodes × 40 tasks = 120 ranks.
- **Prior working job:** `moea_09_reeval` 208863, **same code**, 16 ranks — printed `Rank 0: Scattering realization data to 16 ranks...` and proceeded past `load()` without incident.
- **Logs (full context):**
  - `/home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/MOEA-FIND/workflows/slurm/slurm_logs/moea_09_reeval_209166.err` (120 ranks — failure)
  - `/home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/MOEA-FIND/workflows/slurm/slurm_logs/moea_09_reeval_209166.out`
  - `/home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/MOEA-FIND/workflows/slurm/slurm_logs/moea_09_reeval_208863.err` (16 ranks — past `load()`)

**Root cause.** The current implementation has rank 0 read every realization from the ensemble HDF5 into a giant dict (`all_data`), then pickle-scatter per-rank slices to every other rank. With 1446 realizations × ~20 years daily × dozens of nodes, the aggregate pickled payload rank 0 tries to emit exceeds what mpi4py's pickle-based `scatter` can reliably serialize — specifically, either the aggregate send buffer trips an `INT_MAX` byte-count check, or the OpenMPI eager-send path rejects the combined payload. The `MPI_ERR_ARG` comes from the C-level argument validation in `PMPI_Scatter`. This is an inherent limitation of pickle-based collectives at high rank counts, not an OpenMPI/mpi4py bug we can silence.

`np.array_split`ing across 16 ranks keeps each slice and the aggregate payload small enough to fit; 120 ranks does not. The fix must remove the scatter, not tune it.

---

## 2. File and code location

File: `src/pywrdrb/pre/predict_inflows.py`
Class: `PredictedInflowEnsemblePreprocessor` (inherits `PredictedInflowPreprocessor`)
Method to modify: `load()` — currently lines **479–548**.
Related method (verify it still works): `process()` — lines 549–598. It consumes `self.realization_data` via `self.realization_data[str(realization_id)]` for the rank's own slice, so any fix must leave `self.realization_data` populated with the current rank's realizations, keyed by string id.

Current `load()` for reference (as of the handoff commit):

```python
def load(self):
    """Load available realization IDs, catchment water consumption, and all realization data.

    Rank 0 performs all file reads. Data is distributed using MPI primitives
    that avoid pickle-based large-object broadcasts:
    - wc CSV is broadcast as a raw UTF-8 string (avoids DataFrame pickle)
    - realization DataFrames are scattered so each rank receives only its slice
    """

    ### Load water consumption CSV on rank 0; broadcast as raw string to avoid
    # pickle-based DataFrame bcast which fails on some HPC MPI stacks.
    fname = self.input_dirs["sw_avg_wateruse_pywrdrb_catchments_mgd.csv"]

    if self.rank == 0:
        print(f"Rank 0: Loading catchment water consumption data...")
        with open(fname, "r") as f:
            wc_str = f.read()
    else:
        wc_str = None

    if self.use_mpi:
        wc_str = self.comm.bcast(wc_str, root=0)

    wc = pd.read_csv(io.StringIO(wc_str))
    wc.index = wc["node"]
    self.catchment_wc = wc

    ### Rank 0 reads all realization data from HDF5, then scatters each rank's
    # assigned slice. This avoids (a) concurrent file opens and (b) broadcasting
    # a large dict of all DataFrames to every rank.
    if self.rank == 0:
        print(f"Rank 0: Reading all realization data from HDF5...")
        with h5py.File(self.ensemble_hdf5_file, "r") as f:
            if self.realization_ids is None:
                self.realization_ids = [key for key in f.keys()]
            all_data = {
                rid: self._extract_realization_from_open_file(f, rid)
                for rid in self.realization_ids
            }
    else:
        all_data = None

    if self.use_mpi:
        self.realization_ids = self.comm.bcast(self.realization_ids, root=0)
        if self.rank == 0:
            slices = [
                {rid: all_data[rid] for rid in chunk}
                for chunk in np.array_split(self.realization_ids, self.size)
            ]
            print(f"Rank 0: Scattering realization data to {self.size} ranks...")
        else:
            slices = None
        self.realization_data = self.comm.scatter(slices, root=0)
    else:
        self.realization_data = all_data

    # Initialize STARFIT simulator once if perfect_foresight mode is used
    if "perfect_foresight" in self.modes:
        from pywrdrb.pre.generate_presimulated_releases import (
            STARFITOfflineSimulator,
        )

        self._starfit_simulator = STARFITOfflineSimulator(initial_volume_frac=0.8)
        self._starfit_simulator.load_parameters()

    if self.rank == 0:
        print(
            f"Processing {len(self.realization_ids)} realizations across {self.size} processes"
        )
```

---

## 3. Proposed fix

**Design:** each rank opens the ensemble HDF5 in read-only mode and extracts *only its own slice* using `_extract_realization_from_open_file()` (already on the class). No scatter. The water-consumption CSV broadcast stays as-is (small payload, known-good pattern).

### Why this is safe

- HDF5 with default driver supports many concurrent read-only opens on a shared filesystem (Lustre/NFS/GPFS); no special `parallel=True` build is needed. Read-only opens do not contend with each other.
- The Hopper deployment uses Lustre and this pattern is already used elsewhere in pywrdrb (e.g. every non-MPI `FlowEnsemble` parameter opens the same file per-worker).
- The current method's own docstring rationale ("avoids (a) concurrent file opens") is outdated — the failure mode it is trying to avoid (MPI-IO contention on h5py parallel builds) doesn't apply to the default serial-h5py read-only path MOEA-FIND uses.

### Replacement for lines 506–533

```python
### Determine the realization id list. The caller usually passes
# realization_ids explicitly; otherwise rank 0 enumerates from the
# first node group's column_labels (the HDF5 is stored node-first —
# /<node>/attrs['column_labels'] lists realization ids — so
# f.keys() is node names, NOT realization ids).
if self.rank == 0 and self.realization_ids is None:
    with h5py.File(self.ensemble_hdf5_file, "r") as f:
        first_node = pywrdrb_all_nodes[0]
        labels = f[first_node].attrs["column_labels"]
        self.realization_ids = [str(l) for l in labels]

if self.use_mpi:
    self.realization_ids = self.comm.bcast(self.realization_ids, root=0)

### Each rank reads only its own slice, keyed by string realization id.
# Replaces the prior rank-0-reads-all + pickle-scatter pattern, which
# fails with MPI_ERR_ARG at high rank counts because mpi4py's pickle
# scatter trips INT_MAX buffer limits on large nested dicts
# (see moea_09_reeval_209166.err, 120 ranks × 1446 realizations).
if self.use_mpi:
    my_ids = list(np.array_split(self.realization_ids, self.size)[self.rank])
else:
    my_ids = list(self.realization_ids)

if self.rank == 0:
    print(
        f"Rank 0: each of {self.size} ranks reading its own HDF5 slice "
        f"(~{len(my_ids)} realizations/rank) from {self.ensemble_hdf5_file}"
    )

with h5py.File(self.ensemble_hdf5_file, "r") as f:
    self.realization_data = {
        str(rid): self._extract_realization_from_open_file(f, rid)
        for rid in my_ids
    }
```

Important invariants to preserve:

- `self.realization_ids` is set on **every** rank after this block (either from the caller, or broadcast from rank 0).
- `self.realization_data` is keyed by `str(rid)` because `process()` does `self.realization_data[str(realization_id)]` at line 572.
- In the serial path (`use_mpi=False`), `my_ids` is the full list → behavior matches the old `else: self.realization_data = all_data` branch.
- The `self._starfit_simulator` setup block (lines 535–542) stays unchanged, immediately after the new block.
- The final rank-0 status print at 544–547 stays unchanged.

### What to delete

- The `all_data = {...}` dict construction on rank 0 (its contents no longer need to exist on any single rank).
- The `slices = [...]` list and the `self.comm.scatter(slices, root=0)` call.

### Do not change

- The `wc_str` bcast logic at lines 487–504 — it works and is a separate concern.
- The `process()`, `save()`, or `__init__` methods.
- The `_extract_realization_from_open_file` helper (re-used by the new code).

---

## 4. Verification

**Heads-up on the user's environment:** the user is on a shared HPC (Cornell Hopper). Per their standing feedback, do **not** run heavy compute or MPI jobs on the login node — submit via SLURM. Small unit tests (<~1 min, <2 GB) are fine to run locally.

### 4a. Unit test (add to pywrdrb's test suite)

Add a new test file under `tests/` (the repo already has `tests/conftest.py`, `tests/test_data_retrieval.py`, etc.). Suggested: `tests/test_predicted_inflow_ensemble_load.py`. The test should:

1. Build a small synthetic ensemble HDF5 (~4 realizations, 30 days, handful of pywrdrb nodes) in a tmp dir using the same `stored_by_node=True` layout `FlowEnsemble` / `Ensemble.to_hdf5()` produces.
2. Instantiate `PredictedInflowEnsemblePreprocessor(..., use_mpi=False)` and call `load()`. Assert `len(self.realization_data) == 4` and keys are `{"0","1","2","3"}` as strings.
3. Assert `self.realization_ids == ["0","1","2","3"]` (whatever the caller passes round-trips).
4. (Optional) Use `mpi4py` + `MPI.COMM_SELF` or subprocess a 2-rank mpirun via `subprocess.run` to exercise the MPI path, but only if the repo's existing test runner already supports MPI (check `tests/conftest.py`); otherwise leave the MPI path for smoke testing in §4b.

Run locally from the repo root:

```bash
cd /home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/Pywr-DRB
source ../MOEA-FIND/venv/bin/activate   # editable install target
pytest tests/test_predicted_inflow_ensemble_load.py -v
```

### 4b. MPI smoke test (submit via SLURM, not login node)

Use MOEA-FIND's existing dev-ensemble smoke path. From `/home/fs02/pmr82_0001/tja73/Research/DRB/Pywr-DRB/MOEA-FIND`:

```bash
# 2-node, 80-rank smoke on 10 realizations; skips Stage 1 to exercise only
# the MPI path in load()/process(). Dev subset JSON is already on disk.
sbatch --export=ALL,SUBSET=outputs/data_cache/dev_10.json,SKIP_GENERATE=1 \
    workflows/slurm/09_drb_policy_reeval.slurm \
    outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json
```

*Before submitting*, temporarily edit `workflows/slurm/09_drb_policy_reeval.slurm` to set `#SBATCH --nodes=2 --ntasks-per-node=40` (80 ranks) and `--time=00:15:00` so the smoke completes fast and exits at a rank count known to have stressed the old scatter. Pass criterion: logs print `Rank 0: each of 80 ranks reading its own HDF5 slice (~... realizations/rank) from ...` and `prep_predicted_inflows` completes without MPI errors. Revert those slurm edits after the smoke.

The user will run the 120-rank × 1446-realization full verification themselves after the handoff returns — do not submit that yourself.

### 4c. Non-MPI regression

Quickly confirm the serial path still works:

```bash
python -c "
from pywrdrb.pre import PredictedInflowEnsemblePreprocessor
p = PredictedInflowEnsemblePreprocessor(
    flow_type='nhmv10',  # or any registered flow_type with ensemble HDF5
    ensemble_hdf5_file='<path to a small test HDF5>',
    realization_ids=['0','1'],
    use_mpi=False,
)
p.load()
print(sorted(p.realization_data.keys()))
"
```

---

## 5. Commit / PR guidance

- Target branch: whatever the user uses for pywrdrb development (check `git status` / `git log` in the repo; per MOEA-FIND's memory the install is editable from `../Pywr-DRB`, so commit there).
- Commit title: something like `fix(MPI): per-rank HDF5 read in PredictedInflowEnsemblePreprocessor.load()`.
- Commit body should explain: failing configuration (120 ranks × 1446 realizations), the `MPI_ERR_ARG` traceback, and that the fix eliminates the scatter entirely — keeping the non-MPI path identical and bounded to the same wall-clock as before for small rank counts.
- Reference the failing MOEA-FIND job id (`209166`) in the body for traceability.
- Do **not** push to remote without the user's go-ahead (per the user's standing git safety preference — ask before pushing).

---

## 6. Out of scope / do not do

- Do not touch MOEA-FIND. The fix lives entirely in pywrdrb. MOEA-FIND will resubmit against the editable install once the patch lands.
- Do not parallelize HDF5 writes (`save()`) — not in scope and unrelated to the failure.
- Do not modify `PredictedDiversionEnsemblePreprocessor` or `ExtrapolatedDiversionEnsemblePreprocessor` unless they share the same scatter pattern *and* a user check confirms they need the same fix. (Check for `comm.scatter` in their `load()` first; if not present, leave them.)
- Do not rewrite the `wc_str` broadcast or the serial fallback. Minimum viable diff.
- Do not add a `parallel h5py` (MPI-IO) dependency — the fix works with the standard serial h5py build already in the env.

---

## 7. Handback

When finished, reply with:

1. The commit hash and diff summary.
2. Unit test results (pass/fail, any flakiness).
3. Smoke-test job id and a 5-line excerpt from its log showing the new "each of N ranks reading its own HDF5 slice" message and successful Stage 2 completion.
4. Any concerns or follow-ups (e.g. if `process()` turns out to have a related scatter/gather, or if the unit test reveals an edge case).

The user will then resume the MOEA-FIND pipeline restart plan documented at
`/home/fs02/pmr82_0001/tja73/.claude/plans/yesterday-i-was-hopping-logical-moonbeam.md`.
