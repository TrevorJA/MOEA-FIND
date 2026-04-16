#!/bin/bash
#SBATCH --job-name=moea_02_analytic_3d
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Manuscript: §5 Analytic Validation (Figure 2, DD-10)
# Experiment:  02 — 3D analytic proof-of-concept
# Parallelism: serial

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/experiments/02_analytic_3d.py \
    --nfe 50000 \
    --seed 42 \
    --epsilon 0.15 \
    --plot
