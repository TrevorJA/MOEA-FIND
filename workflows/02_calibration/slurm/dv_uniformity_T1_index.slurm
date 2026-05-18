#!/bin/bash
#SBATCH --job-name=02_dv_uniformity_T1_index
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# Stage 02 / DV-uniformity calibration for T=1 INDEX mode (DD-15 A/B test).
# Adds cannonsville_T1_{ad,ks,l2_star} entries to
# outputs/02_calibration/dv_uniformity_calibration/index/calibrated_dv_tolerances.json.
# Index mode at T=1 has n_dvs = 12*(T+1) = 24.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/dv_uniformity_T1_index] calibrating T=1 index mode (n_dvs=24)"
python workflows/02_calibration/dv_uniformity_calibration.py \
    --T 1 \
    --n-boot 2000 \
    --dv-mode index \
    --statistic all \
    --seed 20260420 \
    --site-label cannonsville
