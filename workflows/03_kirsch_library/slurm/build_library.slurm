#!/bin/bash
#SBATCH --job-name=moea_05_kirsch_library
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclusive
# Experiment 05: Generate Kirsch baseline library.
# Serial — SynHydro handles ensemble generation internally.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

: "${LIBRARY_N_TRACES:=10000}"
: "${LIBRARY_N_YEARS:=20}"

python workflows/03_kirsch_library/build_library.py \
    --n-traces "${LIBRARY_N_TRACES}" \
    --n-years "${LIBRARY_N_YEARS}" \
    --ssi 3 \
    --seed 42

echo "exit: $?"
