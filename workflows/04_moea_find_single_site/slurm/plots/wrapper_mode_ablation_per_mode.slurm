#!/bin/bash
#SBATCH --job-name=04_plots_wm_mode
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G
#SBATCH --array=0-5
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup
MODES=(index index index residual residual residual)
SLUGS=(
  index_T20_nfe200000_s42_constrained
  index_T20_nfe200000_s101_constrained
  index_T20_nfe200000_s2026_constrained
  residual_T20_nfe200000_s42_constrained
  residual_T20_nfe200000_s101_constrained
  residual_T20_nfe200000_s2026_constrained
)
IDX="${SLURM_ARRAY_TASK_ID}"
python src/plotting/04_moea_find_single_site/wrapper_mode_ablation.py \
    --mode "${MODES[$IDX]}" --slug "${SLUGS[$IDX]}"
