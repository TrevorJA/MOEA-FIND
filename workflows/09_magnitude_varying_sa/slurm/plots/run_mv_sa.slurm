#!/bin/bash
#SBATCH --job-name=09_plots_run_mv_sa
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --account=tja73
#SBATCH --partition=normal
# Stage 09 / plots/run_mv_sa -- render MV-SA figures from cached parquet.
# Pair with: workflows/09_magnitude_varying_sa/slurm/run_mv_sa.slurm
#
# Slug must match the compute driver. Edit the SLUG variable to point at
# the variant whose figures you want to (re)render.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

SLUG="mvsa__src=residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad__axis=nyc_min_storage_frac__resp=within_trace_percentile__metric_set=h30e96e__methods=delta__n_perc=19__s=42"

python src/plotting/09_magnitude_varying_sa/run_mv_sa.py --slug "${SLUG}"
