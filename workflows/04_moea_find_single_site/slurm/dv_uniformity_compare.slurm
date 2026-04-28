#!/bin/bash
#SBATCH --job-name=moea_14_dv_uniformity_compare
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# DV-uniformity ablation: post-hoc comparison and SI figure rendering.
# Depends on 13_dv_uniformity_ablation.slurm output under
# outputs/exp13_dv_uniformity_ablation/{hydrologic,dv_uniform}/.
#
# Typical chain:
#   ABLATION_JOBID=$(sbatch --parsable --export=ALL,BORG_NFE=500000 \
#                    --array=0-5 workflows/04_moea_find_single_site/slurm/dv_uniformity_ablation.slurm)
#   sbatch --dependency=afterok:${ABLATION_JOBID} \
#          workflows/04_moea_find_single_site/slurm/dv_uniformity_compare.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

python workflows/04_moea_find_single_site/dv_uniformity_compare.py
