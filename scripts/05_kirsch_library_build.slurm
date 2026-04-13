#!/bin/bash
#SBATCH --job-name=moea_05_kirsch_library_build
#SBATCH --output=scripts/slurm_logs/%x_%j.out
#SBATCH --error=scripts/slurm_logs/%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
# Manuscript: §6.3 Library Baseline (Figure 7 baseline)
# Experiment:  05 — Kirsch library generation
# Parallelism: mpi4py work queue; trace IDs farmed round-robin across ranks

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

: "${LIBRARY_N_TRACES:=10000}"
: "${LIBRARY_N_YEARS:=15}"

moea_mpi_launch "${SLURM_NTASKS}" python scripts/05_kirsch_library_build.py \
    --n-traces "${LIBRARY_N_TRACES}" \
    --n-years "${LIBRARY_N_YEARS}" \
    --ssi 3 \
    --seed 42
