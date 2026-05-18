#!/bin/bash
#SBATCH --job-name=02_t_sensitivity_historical
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-3
# Stage 02 / t_sensitivity_historical -- per-T historical-block screening.
# Four-task array indexes T = {5, 10, 20, 30} (coarse-bracket grid).
# Each task writes outputs/02_calibration/t_sensitivity_historical/T{T:02d}/.
# Run t_sensitivity_aggregate.slurm after all tasks complete.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

T_GRID=(5 10 20 30)
T_YEARS="${T_GRID[${SLURM_ARRAY_TASK_ID}]}"

echo "[02/t_sensitivity_historical] T_years=${T_YEARS}"
python workflows/02_calibration/t_sensitivity_historical.py \
    --T-years "${T_YEARS}" \
    --stride 1 \
    --correlation-cap 0.6
