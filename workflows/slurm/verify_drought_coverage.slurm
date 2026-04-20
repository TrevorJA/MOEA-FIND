#!/bin/bash
#SBATCH --job-name=moea_verify_drought_coverage
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Stage A — pass/fail verification (Criteria 1–5) against a completed
# MOEA-FIND archive. Serial; one process reads results.json, applies the
# criteria, writes verification_report.json + five figures.
#
# Usage:
#   sbatch workflows/slurm/verify_drought_coverage.slurm <path/to/results.json>
#
# Example (200K NFE Cannonsville archive):
#   sbatch workflows/slurm/verify_drought_coverage.slurm \
#       outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

RESULTS_JSON="${1:?Usage: sbatch verify_drought_coverage.slurm <path/to/results.json>}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python workflows/diagnostics/verify_drought_coverage.py "${RESULTS_JSON}"
