#!/bin/bash
#SBATCH --job-name=moea_diag_dv_uniformity_calibration
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Diagnostic: bootstrap-calibrate DV-space uniformity tolerance (ablation SI).
# Serial — no Kirsch generator calls, only numpy U[0,1] draws and scalar stats.
# Must complete before workflows/slurm/13_dv_uniformity_ablation.slurm.
#
# Override T, N_BOOT, DV_MODE, STATISTIC, or SEED via the sbatch command line:
#   sbatch --export=ALL,T=20,N_BOOT=2000 \
#          workflows/slurm/diag_dv_uniformity_calibration.slurm

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

: "${T:=20}"
: "${N_BOOT:=2000}"
: "${DV_MODE:=residual}"
: "${STATISTIC:=both}"
: "${SEED:=20260420}"

export OMP_NUM_THREADS=1

python workflows/diagnostics/diag_dv_uniformity_calibration.py \
    --T "${T}" \
    --n-boot "${N_BOOT}" \
    --dv-mode "${DV_MODE}" \
    --statistic "${STATISTIC}" \
    --seed "${SEED}"
