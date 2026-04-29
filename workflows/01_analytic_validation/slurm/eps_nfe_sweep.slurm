#!/bin/bash
#SBATCH --job-name=01_eps_nfe_sweep
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-35%8
# Stage 01 / eps_nfe_sweep -- compute only (MM Borg via MPI).
# Manuscript: Fig 3.
# Grid = 4 epsilons x 3 NFEs x 3 seeds = 36 cells. Each cell = 16 MPI ranks;
# borg_runner auto-picks n_islands=2 (1 controller + 2 islands x (1 master +
# 6 workers)). At n=8 ranks some hard cells (eps=0.10 nfe=50000 seed=91) still
# flake; at n=16 they pass deterministically. Throttled to 8 concurrent cells
# (%8) so 36 mpiruns don't fire simultaneously and trigger Borg's per-node
# startup race. Concurrent peak 8 x 16 = 128 ranks.
# (Borg MM with n_islands=1 at >=7 workers also segfaults; borg_runner._auto_islands
# enforces >=2 islands at n>=5.)
# After all array tasks finish, aggregate via:
#     python workflows/01_analytic_validation/eps_nfe_sweep.py --mode aggregate
# Then plot via:
#     sbatch workflows/01_analytic_validation/slurm/plots/eps_nfe_sweep.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

moea_mpi_launch "${SLURM_NTASKS}" python workflows/01_analytic_validation/eps_nfe_sweep.py \
    --mode cell \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --metric-cap 0
