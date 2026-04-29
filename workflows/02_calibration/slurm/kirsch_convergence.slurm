#!/bin/bash
#SBATCH --job-name=02_kirsch_convergence
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / kirsch_convergence -- compute only.
# SI-D. Reads stage 03 library; optionally a stage 04 Pareto archive.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/02_calibration/kirsch_convergence.py \
    --library outputs/03_kirsch_library/build_library/n10000_t20_ssi3_s42/characteristics.npz \
    --moea-front outputs/04_moea_find_single_site/run_moea_find/residual_T20_nfe200000_s42_constrained/results.json \
    --n-seeds 10 \
    --seed 42
