#!/bin/bash
#SBATCH --job-name=03_recharacterize_library_T5
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
# DD-15 follow-up: re-extract metrics from the existing T=5 Kirsch library
# (no regeneration) so the Tier-G hazard-clean metrics land in the
# `characteristics_extended.npz` for Stage-2 KS comparison and Stage-3
# decision matrix consumption.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/03_kirsch_library/recharacterize_library.py \
    --T-years 5 \
    --n-traces 10000 \
    --seed 42
