#!/bin/bash
#SBATCH --job-name=02_plots_dv_uniformity_calibration
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-1
# Stage 02 / plots / dv_uniformity_calibration -- two tasks (residual + index).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

MODES=(residual index)
DV_MODE="${MODES[${SLURM_ARRAY_TASK_ID}]}"

python src/plotting/02_calibration/dv_uniformity_calibration.py --dv-mode "${DV_MODE}"
