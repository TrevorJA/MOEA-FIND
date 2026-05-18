#!/bin/bash
#SBATCH --job-name=01_plots_analytic_2d
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Stage 01 / plots / analytic_2d -- read cached pareto.npz, render PDF.
# Slug must match the compute-driver invocation in
# workflows/01_analytic_validation/slurm/analytic_2d.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/01_analytic_validation/analytic_2d.py \
    --slug k2_nfe100000_eps0.060_s42 \
    --seed 42
