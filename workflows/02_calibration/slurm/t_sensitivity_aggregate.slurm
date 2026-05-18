#!/bin/bash
#SBATCH --job-name=02_t_sensitivity_aggregate
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / t_sensitivity_aggregate -- cross-T aggregation. Run AFTER
# t_sensitivity_historical.slurm completes for the full T-grid.
# Submit with:
#   sbatch --dependency=afterok:<historical_job_id> \
#          workflows/02_calibration/slurm/t_sensitivity_aggregate.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/t_sensitivity_aggregate] aggregating T-grid"
python workflows/02_calibration/t_sensitivity_aggregate.py \
    --T-grid 5 10 20 30
