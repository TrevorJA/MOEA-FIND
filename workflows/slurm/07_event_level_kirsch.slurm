#!/bin/bash
#SBATCH --job-name=moea_07_event_level_kirsch
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
# Manuscript: §6.4 Event-level Formulation
# Experiment:  07 — short-trace event-level Kirsch MOEA-FIND
# Parallelism: MPI master-worker (MM Borg)
# STATUS:      SCAFFOLDED; enable once event-level helpers are merged.

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

moea_mpi_launch "${SLURM_NTASKS}" python workflows/experiments/07_event_level_kirsch.py \
    --nfe 50000 \
    --n-years 8 \
    --ssi 3 \
    --seed 42
