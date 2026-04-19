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

# Resolve which exp04 variant feeds fig05. If BORG_NFE is set (e.g. by
# the upstream exp04 SLURM via --export), pass it through; otherwise
# let script 10 fall back to DEFAULT_EXPERIMENT.nfe.
PY_ARGS=()
if [ -n "${BORG_NFE:-}" ]; then
    PY_ARGS+=(--exp04-nfe "${BORG_NFE}")
fi

python workflows/experiments/10_plot_manuscript_figures.py "${PY_ARGS[@]}"
