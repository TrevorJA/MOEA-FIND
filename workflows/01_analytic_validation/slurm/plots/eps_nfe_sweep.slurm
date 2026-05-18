#!/bin/bash
#SBATCH --job-name=01_plots_eps_nfe_sweep
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Stage 01 / plots / eps_nfe_sweep -- read aggregate.json, render heatmap.
# Make sure aggregate.json exists first:
#     python workflows/01_analytic_validation/eps_nfe_sweep.py --mode aggregate

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/01_analytic_validation/eps_nfe_sweep.py
