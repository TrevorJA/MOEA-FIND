#!/bin/bash
#SBATCH --job-name=08_plots_compare_methods
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Stage 08 / plots / compare_methods -- side-by-side SA across run_sa slugs.
# Edit --slugs / --labels / --method / --tag below to match the comparison.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python src/plotting/08_nyc_sensitivity/compare_methods.py \
    --slugs sa__src=residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad__metric_set=h30e96e__methods=delta-pawn-rbd_fast__n_outcomes=4__s=42 \
    --labels primary \
    --method delta \
    --tag default
