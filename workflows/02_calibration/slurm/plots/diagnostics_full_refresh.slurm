#!/bin/bash
#SBATCH --job-name=02_diagnostics_full_refresh
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Refresh DF1, DF2, DF3, DF4, DF5, DF6 against the latest decision_matrix
python src/plotting/02_calibration/diagnostics_for_decision.py
python src/plotting/02_calibration/diagnostics_hazard_clean_corr.py
python src/plotting/02_calibration/diagnostics_corr_compare.py
