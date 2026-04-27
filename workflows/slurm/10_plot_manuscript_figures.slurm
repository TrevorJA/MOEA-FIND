#!/bin/bash
#SBATCH --job-name=moea_10_manuscript_figures
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
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
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Resolve which exp04 variant feeds fig05/06. If BORG_NFE is set
# (e.g. by the upstream exp04 SLURM via --export), pass it through;
# otherwise default to the production NFE used by current exp04/09
# runs.  DEFAULT_EXPERIMENT.nfe in src/experiment_config.py is a
# plumbing smoke value (nfe=100) that would cause fig05/06 to SKIP
# with "missing outputs/exp04_.../residual_T20_nfe100_s42_constrained".
: "${BORG_NFE:=200000}"
python workflows/experiments/10_plot_manuscript_figures.py \
    --exp04-nfe "${BORG_NFE}"
