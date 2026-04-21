#!/bin/bash
#SBATCH --job-name=moea_diag_kirsch_wrapper_fidelity
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
# Diagnostic: Kirsch-wrapper fidelity under uniform random DVs + drought-space
# coverage. Compares residual vs index vs kirsch_baseline vs historical at the
# trace level (Phase A) and in drought-characteristic space (Phase B).
#
# Phase A answers whether each wrapper reproduces the underlying Kirsch
# distribution when DVs are drawn uniformly. Phase B overlays mode clouds on
# the exp05 10k Kirsch library as the extended reference. Required to finish
# cleanly before Phase D (scripts 15/16) is considered.
#
# Overrides (via --export=ALL,...):
#   N_TRACES  (default 1000)   traces per ensemble
#   N_YEARS   (default 20)     trace length in water years
#   SSI       (default 3)      SSI accumulation period
#   SEED      (default 42)     random seed for baseline + wrapper DV draws
#
# Usage:
#   sbatch workflows/slurm/diag_kirsch_wrapper_fidelity.slurm
#   sbatch --export=ALL,N_TRACES=2000,SEED=101 \
#          workflows/slurm/diag_kirsch_wrapper_fidelity.slurm

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

: "${N_TRACES:=1000}"
: "${N_YEARS:=20}"
: "${SSI:=3}"
: "${SEED:=42}"

export OMP_NUM_THREADS=1

python workflows/diagnostics/diag_kirsch_wrapper_fidelity.py \
    --n-traces "${N_TRACES}" \
    --n-years "${N_YEARS}" \
    --ssi "${SSI}" \
    --seed "${SEED}"
