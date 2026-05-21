#!/bin/bash
#SBATCH --job-name=08_run_sa
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --account=tja73
#SBATCH --partition=normal
# Stage 08 / run_sa -- compute only.
# Manuscript: SI-12 (NYC sensitivity to drought-hazard characteristics).
#
# Bank/chars paths are baked: change them by editing this file, not via --export.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

# Production first-event slug. Override by setting MOEA_SLUG, or pass
# explicit --bank/--chars via $1/$2 if you need different paths.
MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
BANK="${1:-outputs/06_pywrdrb_reeval/policy_reeval/${MOEA_SLUG}/results/metric_bank.parquet}"
CHARS="${2:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

python workflows/08_nyc_sensitivity/run_sa.py \
    --bank  "${BANK}" \
    --chars "${CHARS}" \
    --config workflows/08_nyc_sensitivity/configs/delta.yaml \
    --methods delta \
    --seed 42 \
    --n-bootstrap 10 \
    --n-convergence-replicates 0 \
    --convergence-sizes 1000 \
    --delta-num-resamples 10
