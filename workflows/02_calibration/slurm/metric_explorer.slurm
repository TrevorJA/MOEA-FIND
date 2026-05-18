#!/bin/bash
#SBATCH --job-name=02_metric_explorer
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / metric_explorer -- exploratory drought-metric selection.
# Computes ~28-metric candidate library on historical T-blocks, runs
# spread / correlation / clustering / PCA, writes recommendation.md
# plus diagnostic figures.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/metric_explorer.py \
    --T-years 20 \
    --stride 1 \
    --correlation-cap 0.6

python src/plotting/02_calibration/metric_explorer.py
