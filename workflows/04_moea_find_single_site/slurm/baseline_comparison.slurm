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

python workflows/04_moea_find_single_site/baseline_comparison.py \
    --kirsch-library outputs/03_kirsch_library/build_library/n10000_t20_ssi3_s42/characteristics.npz \
    --moea-front outputs/04_moea_find_single_site/run_moea_find/residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad/results.json
