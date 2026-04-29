#!/bin/bash
#SBATCH --job-name=02_plots_metric_blocks
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 02 / plots / metric_blocks -- read cached compute outputs, render PDFs.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/02_calibration/plots/metric_blocks.py
