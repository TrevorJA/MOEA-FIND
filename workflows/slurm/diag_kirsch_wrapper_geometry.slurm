#!/bin/bash
#SBATCH --job-name=moea_diag_kirsch_wrapper_geometry
#SBATCH --output=workflows/slurm/slurm_logs/%x_%A.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%A.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8G
# Diagnostic: Kirsch-wrapper mapping geometry (index vs residual).
# Sweeps a single DV coordinate over [0,1] for both wrapper modes and
# produces two SI figures: fig_mapping_smoothness.pdf, fig_mapping_cartoon.pdf.
#
# Override SEED, N_SWEEP, or N_YEARS via the sbatch command line:
#   sbatch --export=ALL,SEED=0,N_SWEEP=800 \
#          workflows/slurm/diag_kirsch_wrapper_geometry.slurm

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

: "${SEED:=42}"
: "${N_SWEEP:=400}"
: "${N_YEARS:=20}"

export OMP_NUM_THREADS=1

python workflows/diagnostics/diag_kirsch_wrapper_geometry.py \
    --seed "${SEED}" \
    --n-sweep-points "${N_SWEEP}" \
    --n-years "${N_YEARS}"
