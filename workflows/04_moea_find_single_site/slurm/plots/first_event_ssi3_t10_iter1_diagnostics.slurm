#!/bin/bash
#SBATCH --job-name=04_plots_first_event_ssi3_t10_iter1_diagnostics
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=0-3
# Re-run the existing K-agnostic archive_diagnostics figure on every
# iter1 seed of the first-event SSI-3 family. Pairwise Pareto scatter,
# trace montage, DV-uniformity AD constraint distribution, hyperplane
# sum check, per-month flow CDF vs historical. Submit with
# --dependency=afterok:<iter1_array_jobid>.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

SEEDS=(42 7 13 19)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
SEED="${SEEDS[$IDX]}"
SLUG="residual_T10_nfe200000_s${SEED}_constrained_cmdv_uniform_sfxiter1_stad"

echo "[plots/04/iter1_diagnostics] array=${IDX} seed=${SEED} slug=${SLUG}"
python src/plotting/04_moea_find_single_site/archive_diagnostics.py \
    --slug "${SLUG}"
