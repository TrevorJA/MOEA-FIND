#!/bin/bash
#SBATCH --job-name=06_aggregate_windowed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
# Stage 06 (WINDOWED) aggregation -- run ONCE after every
# policy_reeval_windowed array task COMPLETES. Concatenates per-scenario
# outputs into results/metric_bank.parquet (+ .csv), satisficing table,
# and window_manifest.parquet for stages 07/08/09. Serial, single node.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
PARETO_RESULTS="${1:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

export PYTHONUNBUFFERED=1
echo "[06w-agg.slurm] Start: $(date)"

python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/aggregate_windowed.py" \
    --pareto-results "${PARETO_RESULTS}"

echo "[06w-agg.slurm] End: $(date)"
