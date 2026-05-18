#!/bin/bash
#SBATCH --job-name=06_policy_reeval_prep_fixed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal

# Stage 06 / PREP ONLY -- single-node isolation of the pywrdrb preprocessors.
#
# Why a separate 1-node job: jobs 229834/229894 hung and 229956/229963
# failed because Hopper's inter-node IB/TCP fabric is unreliable. Rank-
# placement tricks (-np 20 inside a 10-node allocation) were insufficient:
# OpenMPI's ORTE still launches an orted daemon on every allocated node and
# the mpirun<->orted OOB channel runs over the flaky network, producing
# "unable to complete a TCP connection" and an empty presimulated_releases
# HDF5 (all reservoirs missing). A nodes=1 allocation removes ALL inter-node
# fabric dependency: ORTE spans one node, MPI uses shared-memory vader only.
#
# Generate (Stage 1) is SKIPPED -- its outputs persist for this NEW slug
# (gage_flow_mgd.hdf5, catchment_inflow_mgd.hdf5). This job runs only
# Stage 2 (prep): Step 0 STARFITReleaseEnsemblePreprocessor + Step 1
# PredictedInflowEnsemblePreprocessor (bridge prep_predicted_inflows).
# Simulate (Stage 3, needs the memory-tuned 10-node/200-rank topology) and
# classify run in the separate policy_reeval_fixed.slurm afterward.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxfixed_stad/results.json"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

N_PROCS="${SLURM_NTASKS:-20}"

echo "[06.prep.slurm] Pareto results: ${PARETO_RESULTS}"
echo "[06.prep.slurm] Single-node prep, ${N_PROCS} ranks (vader only)"
echo "[06.prep.slurm] Start: $(date)"

# Single node => all ranks intra-node => shared-memory transport only.
# self,vader (no tcp): any stray inter-node attempt fails fast, never hangs.
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca pml ob1 \
    --mca btl self,vader \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage prep \
        --n-years 10 \
        --demand-source constant_max \
        --use-mpi

echo "[06.prep.slurm] End: $(date)"
