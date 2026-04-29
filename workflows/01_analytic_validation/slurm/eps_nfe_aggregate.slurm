#!/bin/bash
#SBATCH --job-name=01_eps_nfe_aggregate
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup
python workflows/01_analytic_validation/eps_nfe_sweep.py --mode aggregate
