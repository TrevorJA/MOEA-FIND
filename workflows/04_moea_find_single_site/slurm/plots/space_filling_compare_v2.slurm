#!/bin/bash
#SBATCH --job-name=04_plots_space_filling_compare_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Space-filling comparison: historical (n=73), random-DV Kirsch (n=10000),
# and Pareto (n=95k+) on the v2 bounded [0,1]^K hypercube. Reports four
# discrepancies + fill distance + min pairwise distance, full-set and
# bootstrap-matched-to-73.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

python src/plotting/04_moea_find_single_site/space_filling_compare.py \
    --slug residual_T1_nfe200000_s42_constrained_cmdv_uniform_sfxv2_stad \
    --n-kirsch 10000 \
    --n-boot 500 \
    --n-grid-fd 4096 \
    --seed 20260506
