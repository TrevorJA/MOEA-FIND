#!/bin/bash
#SBATCH --job-name=02_metric_blocks
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / metric_blocks -- compute only.
# SI-E. Sliding T-year window over the historical record; per-metric variability.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/metric_blocks.py \
    --T-years 20 \
    --stride 1 \
    --ssi 3
