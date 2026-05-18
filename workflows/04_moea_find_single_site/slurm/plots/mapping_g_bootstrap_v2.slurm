#!/bin/bash
#SBATCH --job-name=04_plots_mapping_g_bootstrap_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Mapping G bootstrap-uncertainty SI diagnostic for the v2 archive.
# 2000 bootstrap replicates of (μ̂*, σ̂*) on the n=73 historical record;
# propagates into D_g(s_j) bands + SE(D_g) profile per window.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python src/plotting/04_moea_find_single_site/mapping_g_bootstrap.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad \
    --n-boot 2000 \
    --n-grid 200 \
    --seed 20260506
