#!/bin/bash
#SBATCH --job-name=06_verify_drought_coverage_fixed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Stage 06 / verify_drought_coverage -- pre-flight pass/fail report
# (Criteria 1-5) on a stage 04 archive. Pair with
# workflows/06_pywrdrb_reeval/slurm/plots/verify_drought_coverage.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Bake the Pareto-archive path. To verify a different stage 04 run,
# edit this line; do NOT pass arguments via sbatch.
PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxfixed_stad/results.json"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python -u workflows/06_pywrdrb_reeval/verify_drought_coverage.py \
    --pareto-results "${PARETO_RESULTS}"
