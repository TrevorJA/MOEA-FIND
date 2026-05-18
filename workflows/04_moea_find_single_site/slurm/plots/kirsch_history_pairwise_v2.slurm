#!/bin/bash
#SBATCH --job-name=04_plots_kirsch_history_pairwise_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Companion to bounded_archive_extras: pairwise scatter of N=10000 random-DV
# Kirsch synthetic years vs historical, in the v2 bounded metric space.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python src/plotting/04_moea_find_single_site/kirsch_history_pairwise.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad \
    --n-samples 10000 --seed 20260506
