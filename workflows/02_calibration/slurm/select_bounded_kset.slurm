#!/bin/bash
#SBATCH --job-name=02_select_bounded_kset
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / select_bounded_kset -- DD-15c Phase 3.
# Diagnostic-driven (mapping × K=4 windows) selection over the 48-
# candidate bounded pool. Depends on recompute_bounded_candidates.
# Submit with:
#   sbatch --dependency=afterok:<recompute_jobid> \
#          workflows/02_calibration/slurm/select_bounded_kset.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/select_bounded_kset] running"
python workflows/02_calibration/select_bounded_kset.py
