#!/bin/bash
#SBATCH --job-name=04_baseline_comparison
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --exclusive
# Stage 04 / baseline_comparison -- compute only.
# Manuscript: Fig 7. Reads stage 03 Kirsch library + stage 04 production
# Pareto front; writes comparison_summary.json + pooled.npz.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Production first-event run. The MOEA slug is the one run_moea_find.py
# emits via moea_slug(mode=residual, T=10, nfe=200k, ssi=3,
# metrics=first_event_ssi3_t10, cons=dv-l2, extra={st:ad}, seed=42).
MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
KIRSCH_LIBRARY="${1:-outputs/03_kirsch_library/build_library/library__N=10000__T=20__s=42/characteristics.npz}"
MOEA_FRONT="${2:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

python workflows/04_moea_find_single_site/baseline_comparison.py \
    --kirsch-library "${KIRSCH_LIBRARY}" \
    --moea-front "${MOEA_FRONT}"
