#!/bin/bash
#SBATCH --job-name=02_t_sensitivity_historical_T5
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# DD-15 follow-up: re-run Stage 1 only at T=5 with the Tier-G hazard-clean
# metrics added to src/metrics/extended.py (2026-04-29).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/t_sensitivity_historical.py \
    --T-years 5 \
    --stride 1 \
    --correlation-cap 0.6
