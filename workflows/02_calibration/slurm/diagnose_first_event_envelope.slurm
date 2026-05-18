#!/bin/bash
#SBATCH --job-name=02_diagnose_first_event_envelope
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Phase B + C diagnostic for the T=10y first-SSI3-drought pivot:
#   - Per-event characteristics of every critical SSI-3 drought in the
#     historical record (the 5-axis envelope the MOEA targets).
#   - Random-DV coverage probe at T=10y: empty-event rate (constraint
#     binding rate) + per-month Anderson-Darling plausibility statistics
#     (synthetic vs historical monthly flows).
# Outputs land under outputs/02_calibration/first_event_envelope/ and
# figures/first_event_envelope/.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "[02/diagnose_first_event] T=10, residual, N=2000"
python workflows/02_calibration/diagnose_first_event_envelope.py \
    --n-years 10 \
    --n-samples 2000 \
    --dv-mode residual \
    --ssi-timescale 3 \
    --seed 20260506
