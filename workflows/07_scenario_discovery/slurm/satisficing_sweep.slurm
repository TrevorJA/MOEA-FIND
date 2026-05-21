#!/bin/bash
#SBATCH --job-name=07_satisficing_sweep
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

MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
PARETO_RESULTS="${1:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

python workflows/07_scenario_discovery/satisficing_sweep.py \
    --pareto-results "${PARETO_RESULTS}" \
    --manifest workflows/07_scenario_discovery/satisficing_manifest.yaml \
    --seed 42
