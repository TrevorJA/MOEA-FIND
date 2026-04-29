#!/bin/bash
#SBATCH --job-name=04_plots_dv_uniformity_compare
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/04_moea_find_single_site/plots/dv_uniformity_compare.py
