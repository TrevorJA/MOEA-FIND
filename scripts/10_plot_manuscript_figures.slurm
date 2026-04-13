#!/bin/bash
#SBATCH --job-name=moea_10_manuscript_figures
#SBATCH --output=scripts/slurm_logs/%x_%j.out
#SBATCH --error=scripts/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Manuscript: all figures (§3 concept, §5 validation, §6 single-site,
#             §7 DRB, SI-1 hyperplane, SI-4 convergence)
# Experiment:  10 — final figure assembly
# Parallelism: serial
#
# Run after any experiment output has changed.

# shellcheck source=./_common.sh
source "$(dirname "$0")/_common.sh"
moea_setup

python scripts/10_plot_manuscript_figures.py
