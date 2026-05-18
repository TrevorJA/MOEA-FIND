#!/bin/bash
#SBATCH --job-name=02_dv_uniformity_T1
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Stage 02 / DV-uniformity calibration for T=1 (DD-15).
# Adds cannonsville_T1_{ad,ks,l2_star} entries to
# outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json.
# Must complete before t1_short_block_drb.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/dv_uniformity_T1] calibrating T=1 residual mode (n_dvs=12)"
python workflows/02_calibration/dv_uniformity_calibration.py \
    --T 1 \
    --n-boot 2000 \
    --dv-mode residual \
    --statistic all \
    --seed 20260420 \
    --site-label cannonsville
