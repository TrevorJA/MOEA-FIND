#!/bin/bash
#SBATCH --job-name=04_plots_bounded_archive_extras_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Render the bounded-metric extras (2D Pareto + historical 1-year overlay,
# plus 3D montage) for the short_block_drb_v2 archive at seed 42.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python src/plotting/04_moea_find_single_site/bounded_archive_extras.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad
