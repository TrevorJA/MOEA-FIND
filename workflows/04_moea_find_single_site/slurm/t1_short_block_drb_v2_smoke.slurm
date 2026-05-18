#!/bin/bash
#SBATCH --job-name=04_t1_short_block_drb_v2_smoke
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --exclusive
# Stage 04 / T=1 short-block v2 — smoke wiring test (1k NFE, 4-rank MPI).
# Validates that:
#   * preset 'short_block_drb_v2' resolves through chars_fn to 4 bounded values
#   * anti-ideal D* = 1.0 for all four objectives (CONSTANT rule)
#   * Borg runs to completion and writes results.json
#   * archive_diagnostics is K-agnostic and renders without errors
# Production 200k-NFE run is in t1_short_block_drb_v2.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=${SLURM_NTASKS}
echo "[04/t1_short_block_drb_v2_smoke] tasks=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/t1_short_block_drb_v2.yaml \
    --nfe 1000 \
    --n-islands 1 \
    --slug-suffix smoke
