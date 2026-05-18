#!/bin/bash
#SBATCH --job-name=02_diagnostics_for_decision
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
# Stage 02 / diagnostics_for_decision -- four diagnostic figures
# supporting the (K*, T*) decision under DD-15. Reads Stage 1+2+3
# outputs only; no MOEA archive needed.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/02_calibration/diagnostics_for_decision.py
