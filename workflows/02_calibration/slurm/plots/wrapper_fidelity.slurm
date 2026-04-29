#!/bin/bash
#SBATCH --job-name=02_plots_wrapper_fidelity
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 02 / plots / wrapper_fidelity -- read cached compute outputs, render PDFs.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/02_calibration/plots/wrapper_fidelity.py
