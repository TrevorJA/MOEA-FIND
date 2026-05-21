#!/bin/bash
#SBATCH --job-name=03_build_library
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclusive
# Stage 03 / build_library -- compute only.
# Manuscript: Fig 7 baseline (library generation).
# Serial -- SynHydro handles ensemble generation internally.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/03_kirsch_library/build_library.py \
    --n-traces 10000 \
    --n-years 20 \
    --seed 42
