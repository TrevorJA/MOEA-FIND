#!/bin/bash
#SBATCH --job-name=04_plots_mapping_g_validity_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Mapping G statistical-validity diagnostic for the short_block_drb_v2
# archive (4 windows × 4 SI-quality panels per window). Output PDF/PNG
# and a JSON of test statistics suitable for a SI table.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python src/plotting/04_moea_find_single_site/mapping_g_validity.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad
