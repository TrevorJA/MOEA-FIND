#!/bin/bash
#SBATCH --job-name=04_event_level
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
# Stage 04 / event_level -- compute only.
# SI-H. Short-trace event-level Kirsch MOEA-FIND.
# STATUS: SCAFFOLDED; enable once event-level helpers are merged. Until
# then the driver runs in --dry-run mode.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/04_moea_find_single_site/event_level.py \
    --nfe 50000 \
    --n-years 8 \
    --ssi 3 \
    --seed 42 \
    --dry-run
