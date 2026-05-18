#!/bin/bash
#SBATCH --job-name=07_scenario_discovery_plots_fixed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 07 / scenario_discovery_plots -- final SD figures (main Fig 9).
# Reads stage 06 metric bank + drought levels and stage 04 Pareto archive,
# fits gbt + logreg classifiers, writes single-panel and multi-panel
# satisficing maps plus per-classifier AUC summary bars.
# Serial, no Pywr-DRB, no MPI.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxfixed_stad/results.json"

python src/plotting/07_scenario_discovery/scenario_discovery_plots.py \
    --pareto-results "${PARETO_RESULTS}" \
    --manifest workflows/07_scenario_discovery/satisficing_manifest.yaml \
    --classifier-models gbt logreg \
    --seed 42
