#!/bin/bash
#SBATCH --job-name=moea_diag_constraint_calibration
#SBATCH --output=scripts/slurm_logs/%x_%j.out
#SBATCH --error=scripts/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Diagnostic: bootstrap-calibrate plausibility-constraint tolerances.
# Serial (not MPI) — Kirsch bootstrap + statistic evaluations are vectorised
# enough that one node suffices. Must complete before any run that depends
# on outputs/diag_constraint_calibration/calibrated_tolerances.json.
#
# Override T, N_BOOT, SSI_TIMESCALE, or SEED via the sbatch command line:
#   sbatch --export=ALL,T=20,N_BOOT=2000 scripts/diag_constraint_calibration.slurm

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

: "${T:=20}"
: "${N_BOOT:=2000}"
: "${SSI_TIMESCALE:=3}"
: "${SEED:=20260415}"

export OMP_NUM_THREADS=1

python scripts/diag_constraint_calibration.py \
    --T "${T}" \
    --n-boot "${N_BOOT}" \
    --ssi-timescale "${SSI_TIMESCALE}" \
    --seed "${SEED}"
