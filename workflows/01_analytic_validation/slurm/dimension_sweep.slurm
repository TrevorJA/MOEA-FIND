#!/bin/bash
#SBATCH --job-name=01_dimension_sweep
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --array=2-6
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Stage 01 / dimension_sweep -- compute only (MM Borg via MPI).
# K=2..6 hypercube feasibility. Manuscript: Fig 4.
# Pair with: workflows/01_analytic_validation/slurm/plots/dimension_sweep.slurm
#
# Each cell = 8 MPI ranks; borg_runner auto-picks n_islands=2 (1 controller +
# 2 islands x (1 master + 3 workers)). 5 cells in array, concurrent peak
# ~40 cores -- fits 1 node within the 400-core budget.
# (Borg MM with n_islands=1 at 7 workers segfaults on small allocations;
# borg_runner._auto_islands forces >=2 islands at n>=6.)

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

K="${SLURM_ARRAY_TASK_ID}"

# Per-K epsilon tuned so all K give comparable archive sizes.
case "${K}" in
    2) EPS=0.10 ;;
    3) EPS=0.15 ;;
    4) EPS=0.20 ;;
    5) EPS=0.25 ;;
    6) EPS=0.30 ;;
    *) EPS=0.30 ;;
esac

echo "[01/dimension_sweep] K=${K} eps=${EPS} ntasks=${SLURM_NTASKS}"
moea_mpi_launch "${SLURM_NTASKS}" python workflows/01_analytic_validation/dimension_sweep.py \
    --k "${K}" \
    --nfe 5000 \
    --epsilon "${EPS}" \
    --feasible-shape cube \
    --seed 42
