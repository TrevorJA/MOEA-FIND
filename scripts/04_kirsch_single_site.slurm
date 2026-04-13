#!/bin/bash
#SBATCH --job-name=moea_04_kirsch_single_site
#SBATCH --output=scripts/slurm_logs/%x_%j.out
#SBATCH --error=scripts/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
# Manuscript: §6.1-6.2 Single-Site Kirsch (Figures 5-6)
# Experiment:  04 — single-site Kirsch MOEA-FIND
# Parallelism: MPI master-worker (MM Borg); 1 master + (N-1) evaluators
#
# Override BORG_NFE and BORG_MODE in the sbatch command line if desired:
#     sbatch --export=ALL,BORG_NFE=200000,BORG_MODE=index scripts/04_kirsch_single_site.slurm

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

: "${BORG_NFE:=100000}"
: "${BORG_MODE:=residual}"

moea_mpi_launch "${SLURM_NTASKS}" python scripts/04_kirsch_single_site.py \
    --nfe "${BORG_NFE}" \
    --mode "${BORG_MODE}" \
    --n-years 30 \
    --ssi 3 \
    --seed 42 \
    --plot
