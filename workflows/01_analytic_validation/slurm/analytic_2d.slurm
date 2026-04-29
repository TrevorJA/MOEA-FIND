#!/bin/bash
#SBATCH --job-name=01_analytic_2d
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Stage 01 / analytic_2d -- compute only (MM Borg via MPI).
# Manuscript: Fig 1.
# Pair with: workflows/01_analytic_validation/slurm/plots/analytic_2d.slurm
# Parallelism: 8 MPI ranks; borg_runner picks n_islands=2 on a single node.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

moea_mpi_launch "${SLURM_NTASKS}" python workflows/01_analytic_validation/analytic_2d.py \
    --nfe 100000 \
    --seed 42 \
    --epsilon 0.06
