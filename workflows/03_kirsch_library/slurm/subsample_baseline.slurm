#!/bin/bash
#SBATCH --job-name=moea_06_library_subsample
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Manuscript: §6.3 Library Baseline (Figure 7)
# Experiment:  06 — LHS/Sobol subsample of the Kirsch library
# Parallelism: serial (NumPy/SciPy vectorized)

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

LIBRARY=${LIBRARY:-outputs/exp05_kirsch_library/characteristics.json}

for method in lhs sobol random; do
    python workflows/03_kirsch_library/subsample_baseline.py \
        --library "${LIBRARY}" \
        --n-select 200 \
        --method "${method}" \
        --seed 42
done
