#!/bin/bash
#SBATCH --job-name=02_eval_cost_timing
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-3
# Stage 02 / eval_cost_timing -- per-T MOEA inner-loop wall-time.
# Independent of Stages 1 and 2 (no dependency).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

T_GRID=(5 10 20 30)
T_YEARS="${T_GRID[${SLURM_ARRAY_TASK_ID}]}"

echo "[02/eval_cost_timing] T_years=${T_YEARS}"
python workflows/02_calibration/eval_cost_timing.py \
    --T-years "${T_YEARS}" \
    --mode residual \
    --n-warmup 20 \
    --n-eval 200 \
    --seed 42
