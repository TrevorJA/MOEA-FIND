#!/bin/bash
#SBATCH --job-name=02_wrapper_geometry
#SBATCH --output=workflows/slurm_logs/%x_%A.out
#SBATCH --error=workflows/slurm_logs/%x_%A.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8G
# Stage 02 / wrapper_geometry -- compute only.
# SI-C. Sweeps a single DV coordinate over [0,1] for both wrapper modes;
# produces sweep.npz consumed by the paired plotting driver.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/wrapper_geometry.py \
    --seed 42 \
    --n-sweep-points 400 \
    --n-years 20
