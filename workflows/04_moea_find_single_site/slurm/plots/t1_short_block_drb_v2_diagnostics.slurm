#!/bin/bash
#SBATCH --job-name=04_plots_t1_short_block_drb_v2_diagnostics
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Render the K-agnostic archive_diagnostics PDF for the DD-15c v2 archive.
# Verifies all D_j ∈ [0, 1], no flood-corner Pareto, hyperplane condition holds.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/04_moea_find_single_site/archive_diagnostics.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad
