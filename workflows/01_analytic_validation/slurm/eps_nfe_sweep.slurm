#!/bin/bash
#SBATCH --job-name=moea_03_eps_nfe_sweep
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-89
# Manuscript: §5 Analytic Validation (Figure 3)
# Experiment:  03 — epsilon x NFE sensitivity sweep
# Parallelism: SLURM job array (one task per (eps, nfe, seed) cell)
#
# Default grid = 6 epsilons x 3 NFEs x 5 seeds = 90 cells.
# Update --array if you change DEFAULT_EPSILONS / DEFAULT_NFES / DEFAULT_SEEDS.
#
# Usage:
#     sbatch scripts/03_eps_nfe_sweep.slurm              # run the full array
#     # then, after the array completes:
#     python workflows/01_analytic_validation/eps_nfe_sweep.py --mode aggregate

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/01_analytic_validation/eps_nfe_sweep.py \
    --mode cell \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --metric-cap 0
