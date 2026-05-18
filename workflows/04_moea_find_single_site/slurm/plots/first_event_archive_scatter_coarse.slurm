#!/bin/bash
#SBATCH --job-name=04_fe_archive_scatter_coarse
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Multi-panel 2D scatter matrix + 3D views of the COARSE first-event
# MMBorg archive vs the historical first-event cloud. Light/serial.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup
export PYTHONUNBUFFERED=1

SLUG=residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxcoarse_stad
RES="outputs/04_moea_find_single_site/run_moea_find/${SLUG}/results.json"
OUT="figures/04_moea_find_single_site/archive_scatter/${SLUG}"

echo "[04/fe-scatter] $(date)  results=${RES}"
python src/plotting/04_moea_find_single_site/first_event_archive_scatter.py \
    --results "${RES}" \
    --out-dir "${OUT}" \
    --ssi 3
echo "[04/fe-scatter] end $(date)"
