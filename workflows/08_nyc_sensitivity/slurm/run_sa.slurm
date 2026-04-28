#!/bin/bash
#SBATCH --job-name=moea_08_nyc_sensitivity
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=tja73
#SBATCH --partition=normal
# Manuscript: SI-12 — NYC sensitivity to drought-hazard characteristics
# Experiment: 10 — Stage-08 global SA on Stage-06 NYC outcomes.
# Parallelism: serial (SA index computation is cheap; bootstrap + convergence
# dominate runtime but are still seconds-to-minutes per method).
#
# Inputs (positional, order matters):
#   $1  Stage-06 metric_bank.parquet (Y)
#   $2  Upstream MOEA-FIND results.json with pareto_chars + objective_keys (X)
#
# Optional overrides via --export:
#   SA_CONFIG     YAML preset under workflows/08_nyc_sensitivity/configs/
#                 (default: workflows/08_nyc_sensitivity/configs/all_methods.yaml)
#   SA_SEED       random seed (default: 42)
#   SA_NBOOT      bootstrap replicates (default: 1000)
#
# Usage:
#   sbatch workflows/08_nyc_sensitivity/slurm/run_sa.slurm \
#       outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet \
#       outputs/exp04_kirsch_single_site/<slug>/results.json
#
# Smoke test (small bootstrap + tiny realization subset is configured by
# the YAML config; not by SLURM args):
#   sbatch --time=00:10:00 --export=ALL,SA_NBOOT=50 \
#       workflows/08_nyc_sensitivity/slurm/run_sa.slurm \
#       outputs/exp09_drb_policy_reeval/<slug>/results/metric_bank.parquet \
#       outputs/exp04_kirsch_single_site/<slug>/results.json

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

BANK="${1:?Usage: sbatch run_sa.slurm <metric_bank.parquet> <results.json>}"
CHARS="${2:?Usage: sbatch run_sa.slurm <metric_bank.parquet> <results.json>}"

: "${SA_CONFIG:=workflows/08_nyc_sensitivity/configs/all_methods.yaml}"
: "${SA_SEED:=42}"
: "${SA_NBOOT:=1000}"

PY_ARGS=(
    --bank   "${BANK}"
    --chars  "${CHARS}"
    --config "${SA_CONFIG}"
    --seed   "${SA_SEED}"
    --n-bootstrap "${SA_NBOOT}"
    --plot
)

echo "[08.slurm] bank:   ${BANK}"
echo "[08.slurm] chars:  ${CHARS}"
echo "[08.slurm] config: ${SA_CONFIG}"
echo "[08.slurm] seed=${SA_SEED} n_bootstrap=${SA_NBOOT}"
echo "[08.slurm] start: $(date)"

python "${SLURM_SUBMIT_DIR}/workflows/08_nyc_sensitivity/run_sa.py" \
    "${PY_ARGS[@]}"

echo "[08.slurm] end: $(date)"
