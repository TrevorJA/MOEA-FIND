#!/bin/bash
#SBATCH --job-name=06_plots_verify_drought_coverage
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 06 / plots / verify_drought_coverage -- read cached npz/json,
# render five PDFs. Slug must match the compute-driver invocation in
# workflows/06_pywrdrb_reeval/slurm/verify_drought_coverage.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/06_pywrdrb_reeval/verify_drought_coverage.py \
    --slug residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad
