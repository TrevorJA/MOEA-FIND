#!/bin/bash
#SBATCH --job-name=02_plots_bounded_metric_selection
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / plots / bounded_metric_selection -- DD-15c Phase 4.
# Depends on recompute_bounded_candidates (Phase 2) and
# select_bounded_kset (Phase 3). Submit with:
#   sbatch --dependency=afterok:<select_jobid> \
#          workflows/02_calibration/slurm/plots/bounded_metric_selection.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/plots/bounded_metric_selection] running"
python src/plotting/02_calibration/bounded_metric_selection.py
