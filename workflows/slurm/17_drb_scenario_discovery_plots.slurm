#!/bin/bash
#SBATCH --job-name=moea_17_scenario_discovery
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=tja73
#SBATCH --partition=normal
# Experiment 17 — DRB scenario-discovery figure generator.
# Consumes exp09's persisted metric bank + drought levels and produces
# every scenario-discovery figure (legacy single-panel, multi-panel GBT,
# multi-panel logreg, AUC summary bars).  Serial, no Pywr-DRB, no MPI.
#
# Usage:
#   sbatch workflows/slurm/17_drb_scenario_discovery_plots.slurm \\
#       outputs/exp04_kirsch_single_site/<slug>/results.json
#
# Optional overrides via --export:
#   MANIFEST            satisficing manifest YAML
#                       (default: workflows/experiments/satisficing_manifest.yaml)
#   CLASSIFIERS         space-separated subset of {gbt, logreg}
#                       (default: "gbt logreg")

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

PARETO_RESULTS="${1:?Usage: sbatch 17_drb_scenario_discovery_plots.slurm <path/to/results.json>}"
MODELS="${CLASSIFIERS:-gbt logreg}"

MANIFEST_ARG=()
if [[ -n "${MANIFEST:-}" ]]; then
    MANIFEST_ARG=(--manifest "${MANIFEST}")
fi

echo "[17.slurm] pareto: ${PARETO_RESULTS}"
echo "[17.slurm] models: ${MODELS}"
echo "[17.slurm] start: $(date)"

# shellcheck disable=SC2086
python "${SLURM_SUBMIT_DIR}/workflows/experiments/17_drb_scenario_discovery_plots.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --classifier-models ${MODELS} \
    "${MANIFEST_ARG[@]}"

echo "[17.slurm] end: $(date)"
