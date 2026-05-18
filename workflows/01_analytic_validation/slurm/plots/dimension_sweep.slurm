#!/bin/bash
#SBATCH --job-name=01_plots_dimension_sweep
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 01 / plots / dimension_sweep -- read all k{2..6} dirs, render Fig 4.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/01_analytic_validation/dimension_sweep.py
