#!/bin/bash
#SBATCH --job-name=04_plots_mapping_g_e_diagnostics_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Two SI diagnostics for the v2 archive:
#   1. Mapping E validity (KDE + LOO + body/tail + E-vs-G), 4 panels per window
#   2. Mapping G vs E bootstrap (joint bootstrap, band overlay + SE profile)
# Both are SI-quality and reuse the same archive.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

SLUG=residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad

python src/plotting/04_moea_find_single_site/mapping_e_validity.py --slug "${SLUG}"
python src/plotting/04_moea_find_single_site/mapping_g_e_bootstrap.py \
    --slug "${SLUG}" --n-boot 2000 --n-grid 200 --seed 20260506
