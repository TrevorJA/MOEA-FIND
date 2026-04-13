#!/bin/bash
#SBATCH --job-name=moea_08_drb_multisite
#SBATCH --output=scripts/slurm_logs/%x_%A_%a.out
#SBATCH --error=scripts/slurm_logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-7
# Manuscript: §7.1-7.2 DRB Case Study (Figure 8)
# Experiment:  08 — DRB multi-site MOEA-FIND
# Parallelism: MPI master-worker (MM Borg), one seed per array task
# STATUS:      SCAFFOLDED; enable once multi-site Kirsch wrapper is merged.

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

SEEDS=(42 17 91 5 13 31 77 109)
SEED=${SEEDS[${SLURM_ARRAY_TASK_ID}]}

moea_mpi_launch "${SLURM_NTASKS}" python scripts/08_drb_multisite_moea.py \
    --nfe 500000 \
    --n-years 30 \
    --seeds "${SEED}" \
    --ssi 3
