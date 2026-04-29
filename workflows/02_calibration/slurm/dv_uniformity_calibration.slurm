#!/bin/bash
#SBATCH --job-name=02_dv_uniformity_calibration
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-1
# Stage 02 / dv_uniformity_calibration -- compute only.
# SI-B. Two-task array: index 0 = residual, 1 = index. The driver routes
# outputs to outputs/02_calibration/dv_uniformity_calibration/<mode>/.
# Must complete before stage 04 ablation runs.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

MODES=(residual index)
DV_MODE="${MODES[${SLURM_ARRAY_TASK_ID}]}"

echo "[02/dv_uniformity_calibration] mode=${DV_MODE}"
python workflows/02_calibration/dv_uniformity_calibration.py \
    --T 20 \
    --n-boot 2000 \
    --dv-mode "${DV_MODE}" \
    --statistic both \
    --seed 20260420
