#!/bin/bash
#SBATCH --job-name=moea_09_verify
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal
# One-shot verification of the constant_max preprocessing skip.
# Re-runs Stages 2 + 3 only, on the existing dev_500 subset, so the
# output HDF5 can be compared bit-for-bit against
# /tmp/pywrdrb_output_baseline.hdf5 (snapshotted before the patch).

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

PARETO_RESULTS="outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json"
SUBSET="outputs/data_cache/dev_500.json"
N_PROCS="${SLURM_NTASKS:-16}"

echo "[verify] Pareto: ${PARETO_RESULTS}"
echo "[verify] Subset: ${SUBSET}"
echo "[verify] Processors: ${N_PROCS}"
echo "[verify] Start: $(date)"

# Remove any leftover diversion artifacts from the prior full run so
# the patched simulate step can't accidentally rely on them.
DIV_DIR="outputs/exp09_drb_policy_reeval/residual_T20_nfe200000_s42_constrained/pywrdrb_inputs"
rm -f "${DIV_DIR}/diversion_nj_extrapolated_mgd.hdf5" \
      "${DIV_DIR}/diversion_nyc_extrapolated_mgd.hdf5" \
      "${DIV_DIR}/predicted_diversions_mgd.hdf5" \
      "${DIV_DIR}/predicted_inflows_mgd.hdf5"
echo "[verify] cleaned leftovers from ${DIV_DIR}"
ls -la "${DIV_DIR}"

# Stage 2: prep (MPI)
T2=$(date +%s)
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage prep \
        --n-years 20 \
        --demand-source constant_max \
        --use-mpi \
        --realization-subset-file "${SUBSET}"
T2_END=$(date +%s)
echo "[verify] Stage 2 wall: $((T2_END - T2)) s"

echo "[verify] pywrdrb_inputs/ after Stage 2:"
ls -la "${DIV_DIR}"

# Stage 3: simulate (MPI)
T3=$(date +%s)
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage simulate \
        --n-years 20 \
        --demand-source constant_max \
        --batch-size 10 \
        --use-mpi \
        --realization-subset-file "${SUBSET}"
T3_END=$(date +%s)
echo "[verify] Stage 3 wall: $((T3_END - T3)) s"

echo "[verify] End: $(date)"
