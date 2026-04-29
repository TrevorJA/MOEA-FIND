#!/bin/bash
#SBATCH --job-name=02_wrapper_fidelity
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
# Stage 02 / wrapper_fidelity -- compute only.
# SI-C. Phase A trace-level + Phase B drought-space coverage diagnostics
# for the Kirsch wrapper under uniform random DVs.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/wrapper_fidelity.py \
    --n-traces 1000 \
    --n-years 20 \
    --ssi 3 \
    --seed 42
