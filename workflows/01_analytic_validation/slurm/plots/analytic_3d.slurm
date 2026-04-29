#!/bin/bash
#SBATCH --job-name=01_plots_analytic_3d
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Stage 01 / plots / analytic_3d -- read cached pareto.npz, render PDF.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/01_analytic_validation/plots/analytic_3d.py \
    --slug k3_nfe50000_eps0.150_s42
