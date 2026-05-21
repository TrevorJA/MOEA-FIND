#!/bin/bash
#SBATCH --job-name=06_verify_drought_coverage
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Stage 06 / verify_drought_coverage -- pre-flight pass/fail report
# (Criteria 1-5) on a stage 04 archive. Pair with
# the single workflows/99_supporting_info_figures job (make_si_figures.slurm).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Production first-event slug; pass an alternate results.json path as $1
# or set MOEA_SLUG in the environment to override.
MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
PARETO_RESULTS="${1:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python workflows/06_pywrdrb_reeval/verify_drought_coverage.py \
    --pareto-results "${PARETO_RESULTS}"
