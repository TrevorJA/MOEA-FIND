#!/bin/bash
#SBATCH --job-name=03_subsample_baseline
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Stage 03 / subsample_baseline -- compute only.
# Manuscript: Fig 7 (LHS/Sobol/random library subsamples).
# Loops over the three baseline methods in one job.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

LIBRARY=outputs/03_kirsch_library/build_library/n10000_t20_ssi3_s42/characteristics.json

for method in lhs sobol random; do
    python workflows/03_kirsch_library/subsample_baseline.py \
        --library "${LIBRARY}" \
        --n-select 200 \
        --method "${method}" \
        --seed 42
done
