#!/bin/bash
#SBATCH --job-name=02_constraint_calibration
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / constraint_calibration -- compute only.
# SI-A. Bootstrap-calibrates plausibility-constraint tolerances.
# Must complete before any stage 04 run that consumes
# outputs/02_calibration/constraint_calibration/calibrated_tolerances.json.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/constraint_calibration.py \
    --T 20 \
    --n-boot 2000 \
    --ssi-timescale 3 \
    --seed 20260415
