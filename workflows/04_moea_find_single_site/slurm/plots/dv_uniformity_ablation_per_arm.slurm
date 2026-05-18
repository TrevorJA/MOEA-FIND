#!/bin/bash
#SBATCH --job-name=04_plots_dv_arm
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G
#SBATCH --array=0-7
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup
ARMS=(hydrologic hydrologic hydrologic dv_uniform dv_uniform dv_uniform dv_uniform dv_uniform)
SLUGS=(
  residual_T20_nfe200000_s42_constrained_cmhydrologic
  residual_T20_nfe200000_s101_constrained_cmhydrologic
  residual_T20_nfe200000_s2026_constrained_cmhydrologic
  residual_T20_nfe200000_s42_constrained_cmdv_uniform_stl2_star
  residual_T20_nfe200000_s101_constrained_cmdv_uniform_stl2_star
  residual_T20_nfe200000_s2026_constrained_cmdv_uniform_stl2_star
  residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad
  residual_T20_nfe200000_s101_constrained_cmdv_uniform_stad
)
IDX="${SLURM_ARRAY_TASK_ID}"
python src/plotting/04_moea_find_single_site/dv_uniformity_ablation.py \
    --arm "${ARMS[$IDX]}" --slug "${SLUGS[$IDX]}"
