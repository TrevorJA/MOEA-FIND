#!/bin/bash
#SBATCH --job-name=03_compute_library_metrics
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Stage 03 / compute_library_metrics -- compute only.
# Manuscript: Fig 7 baseline (drought-metric characterization of the
# Kirsch library produced by build_library.slurm).
#
# Runs without regenerating traces; pair with build_library.slurm's
# --n-traces/--n-years/--seed. Pass --ssi to change SSI accumulation
# without touching the library.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/03_kirsch_library/compute_library_metrics.py \
    --n-traces 10000 \
    --n-years 20 \
    --seed 42 \
    --ssi 3
