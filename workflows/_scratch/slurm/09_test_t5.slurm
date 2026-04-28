#!/bin/bash
#SBATCH --job-name=moea_09_t20
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal
# Script 09: Policy re-evaluation with T=20 Pareto front.
# 3 nodes × 40 tasks = 120 MPI ranks for Stages 2-3.
#
# Submit with dependency on Script 04 job:
#   sbatch --dependency=afterok:<SCRIPT04_JOBID> workflows/slurm/09_test_t5.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

PARETO_RESULTS="outputs/exp04_kirsch_single_site/residual_T20_nfe2000000_s42_constrained/results.json"
N_PROCS=$((SLURM_NNODES * SLURM_NTASKS_PER_NODE))

echo "[09] Pareto results: ${PARETO_RESULTS}"
echo "[09] Nodes: ${SLURM_NNODES}, Tasks/node: ${SLURM_NTASKS_PER_NODE}, Total: ${N_PROCS}"
echo "[09] Start: $(date)"

# Verify input exists
if [ ! -f "${PARETO_RESULTS}" ]; then
    echo "[09] ERROR: ${PARETO_RESULTS} not found. Script 04 may not have completed."
    exit 1
fi

# Stage 1 (generate): serial — fits generators, replays DVs
echo "[09] === STAGE 1: GENERATE ==="
python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage generate \
    --n-years 20 \
    --mode residual \
    --seed 42

# Stage 2 (prep): MPI-parallel predicted inflow preprocessing
echo "[09] === STAGE 2: PREP (MPI) ==="
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage prep \
        --n-years 20 \
        --mode residual \
        --seed 42 \
        --use-mpi

# Stage 3 (simulate): MPI-parallel Pywr-DRB batches
echo "[09] === STAGE 3: SIMULATE (MPI) ==="
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage simulate \
        --n-years 20 \
        --mode residual \
        --seed 42 \
        --demand-source constant_max \
        --batch-size 10 \
        --use-mpi

# Stage 4 (classify): serial — extract FFMP levels, scenario discovery
echo "[09] === STAGE 4: CLASSIFY ==="
python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage classify \
    --n-years 20 \
    --mode residual \
    --seed 42 \
    --plot

echo "[09] End: $(date)"
