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
# Pair with: workflows/08_nyc_sensitivity/slurm/plots/run_sa.slurm
#
# Bank/chars paths are baked: change them by editing this file, not via --export.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

python workflows/08_nyc_sensitivity/run_sa.py \
    --bank  outputs/06_pywrdrb_reeval/policy_reeval/residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad/results/metric_bank.parquet \
    --chars outputs/04_moea_find_single_site/run_moea_find/residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad/results.json \
    --config workflows/08_nyc_sensitivity/configs/delta.yaml \
    --methods delta \
    --seed 42 \
    --n-bootstrap 10 \
    --n-convergence-replicates 0 \
    --convergence-sizes 1000 \
    --delta-num-resamples 10
