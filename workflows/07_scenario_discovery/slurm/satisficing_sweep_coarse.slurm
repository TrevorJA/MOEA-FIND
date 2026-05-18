#!/bin/bash
#SBATCH --job-name=07_satisficing_sweep_coarse
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 07 / satisficing_sweep -- cheap GBT manifest sweep.
# Manuscript: SI-J. Consumes stage 04 Pareto archive + stage 06 metric bank.
# No Pywr-DRB, no MPI.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxcoarse_stad/results.json"

python workflows/07_scenario_discovery/satisficing_sweep.py \
    --pareto-results "${PARETO_RESULTS}" \
    --manifest workflows/07_scenario_discovery/satisficing_manifest.yaml \
    --seed 42
