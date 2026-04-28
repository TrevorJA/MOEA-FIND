#!/bin/bash
#SBATCH --job-name=moea_16_wrapper_mode_compare
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
# Wrapper-mode comparison: post-hoc comparison and SI figure rendering.
# Depends on 15_wrapper_mode_ablation.slurm output under
# outputs/exp15_wrapper_mode_ablation/{index,residual}/.
#
# Typical chain:
#   ABLATION_JOBID=$(sbatch --parsable --export=ALL,BORG_NFE=200000 \
#                    --array=0-5 workflows/04_moea_find_single_site/slurm/wrapper_mode_ablation.slurm)
#   sbatch --dependency=afterok:${ABLATION_JOBID} \
#          workflows/04_moea_find_single_site/slurm/wrapper_mode_compare.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/04_moea_find_single_site/wrapper_mode_compare.py
