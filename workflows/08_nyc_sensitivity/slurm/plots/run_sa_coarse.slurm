#!/bin/bash
#SBATCH --job-name=08_plots_run_sa_coarse
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 08 / plots / run_sa -- read cached SA parquets, render diagnostic PDFs.
# Slug must match the compute-driver invocation in
# workflows/08_nyc_sensitivity/slurm/run_sa.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/08_nyc_sensitivity/run_sa.py \
    --slug sa__src=results__metric_set=h7e4695__methods=delta__n_outcomes=4__s=42
