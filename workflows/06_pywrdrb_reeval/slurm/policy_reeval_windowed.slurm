#!/bin/bash
#SBATCH --job-name=06_policy_reeval_windowed
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --array=0-9
# Stage 06 (WINDOWED) per-scenario re-eval of the COARSE stage-04 archive.
#
# Smart hybrid parallelism (per user directive):
#   * SLURM array (10 tasks) chunks the ~1009 scenarios -> robustness +
#     scheduler-friendly; tasks are fully independent.
#   * Within each array task: mpirun -np 30, ONE CORE PER RANK
#     (cpus-per-task=1), each rank owns a disjoint scenario sub-subset
#     and runs the full per-scenario chain serially (batch_size=1).
#   * INTRA-NODE shared-memory only (--mca btl self,vader; nodes=1):
#     ZERO inter-node MPI/ORTE -> structurally immune to the Hopper
#     fabric failures that broke the dense pipeline. No MPI collectives
#     (scenarios independent; aggregation is a separate filesystem pass).
# 10 tasks x 30 ranks = 300 workers over ~1009 scenarios ~= 3-4 each.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxcoarse_stad/results.json"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

N_PROCS="${SLURM_NTASKS:-30}"
echo "[06w.slurm] array task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}, ${N_PROCS} ranks (1 core/rank)"
echo "[06w.slurm] Start: $(date)"

mpirun --oversubscribe -np "${N_PROCS}" \
    --mca pml ob1 \
    --mca btl self,vader \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval_windowed.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --site-label cannonsville \
        --ssi 3 \
        --demand-source constant_max \
        --use-mpi

echo "[06w.slurm] End: $(date)"
