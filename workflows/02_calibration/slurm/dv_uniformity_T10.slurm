#!/bin/bash
#SBATCH --job-name=02_dv_uniformity_T10
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 02 / DV-uniformity calibration for T=10 (first-SSI3-drought pivot).
# Adds cannonsville_T10_{ad,ks,l2_star} entries to
# outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json.
# Required before iter1 production runs of first_event_ssi3_t10.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/dv_uniformity_T10] calibrating T=10 residual mode (n_dvs=120)"
python workflows/02_calibration/dv_uniformity_calibration.py \
    --T 10 \
    --n-boot 2000 \
    --dv-mode residual \
    --statistic all \
    --seed 20260506 \
    --site-label cannonsville
