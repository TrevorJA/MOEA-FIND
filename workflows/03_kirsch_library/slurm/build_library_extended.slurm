#!/bin/bash
#SBATCH --job-name=03_build_library_extended
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=0-3
# Stage 03 / build_library_extended -- 28-metric Kirsch ensemble at every T
# in the joint metric/T sensitivity grid {5, 10, 20, 30}. One task per T;
# each writes outputs/03_kirsch_library/build_library_extended/<slug>/.
# Independent of Stage 1 (no SLURM dependency needed) — generation uses
# only the historical record + seed.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

T_GRID=(5 10 20 30)
T_YEARS="${T_GRID[${SLURM_ARRAY_TASK_ID}]}"

echo "[03/build_library_extended] T_years=${T_YEARS}"
python workflows/03_kirsch_library/build_library_extended.py \
    --n-traces 10000 \
    --n-years "${T_YEARS}" \
    --seed 42
