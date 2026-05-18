#!/bin/bash
#SBATCH --job-name=99_make_figures_fixed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Stage 99 - manuscript figure assembly
# Parallelism: serial
#
# Run after every contributing upstream stage has been executed.

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Single-site MOEA-FIND variant slug that feeds fig05/fig06/fig07.
# Bake the slug here; do not rely on env-var defaults.
SINGLE_SITE_SLUG="residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxfixed_stad"

python workflows/99_manuscript_figures/make_figures.py \
    --single-site-slug "${SINGLE_SITE_SLUG}"
