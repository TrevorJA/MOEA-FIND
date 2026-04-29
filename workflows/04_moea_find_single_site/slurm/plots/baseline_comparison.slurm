#!/bin/bash
#SBATCH --job-name=04_plots_baseline_comparison
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/04_moea_find_single_site/plots/baseline_comparison.py
