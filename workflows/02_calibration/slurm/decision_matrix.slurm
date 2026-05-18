#!/bin/bash
#SBATCH --job-name=02_decision_matrix
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / decision_matrix -- joint K x T decision aggregation.
# Run AFTER all of: t_sensitivity_aggregate, t_sensitivity_kirsch_compare,
# eval_cost_timing array. Submit with:
#   sbatch --dependency=afterok:<agg>:<kirsch_compare>:<timing_array> \
#          workflows/02_calibration/slurm/decision_matrix.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/decision_matrix] running"
python workflows/02_calibration/decision_matrix.py
