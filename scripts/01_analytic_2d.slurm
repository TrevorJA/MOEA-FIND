#!/bin/bash
#SBATCH --job-name=moea_01_analytic_2d
#SBATCH --output=scripts/slurm_logs/%x_%j.out
#SBATCH --error=scripts/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Manuscript: §5 Analytic Validation (Figures 1-2)
# Experiment:  01 — 2D analytic proof-of-concept
# Parallelism: serial

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

python scripts/01_analytic_2d.py \
    --nfe 100000 \
    --seed 42 \
    --plot
