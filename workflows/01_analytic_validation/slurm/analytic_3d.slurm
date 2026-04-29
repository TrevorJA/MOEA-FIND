#!/bin/bash
#SBATCH --job-name=01_analytic_3d
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Stage 01 / analytic_3d -- compute only (MM Borg via MPI).
# Manuscript: Fig 2.
# Pair with: workflows/01_analytic_validation/slurm/plots/analytic_3d.slurm
# Parallelism: 8 MPI ranks; borg_runner picks n_islands=2.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

moea_mpi_launch "${SLURM_NTASKS}" python workflows/01_analytic_validation/analytic_3d.py \
    --nfe 50000 \
    --seed 42 \
    --epsilon 0.15
