#!/bin/bash
#SBATCH --job-name=04_wrapper_mode_ablation
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --array=0-3%3
# Stage 04 / wrapper_mode_ablation -- compute only.
# SI-F. 2 modes x 2 seeds = 4 tasks; throttle %3.
# Array index -> (mode, seed):
#   0: index    s=42
#   1: index    s=101
#   2: residual s=42
#   3: residual s=101
#
# Index mode requires the 252-DV index calibration; residual uses the
# 240-DV residual calibration. Both must exist before submitting.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

MODES=(index index residual residual)
SEEDS=(42 101 42 101)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
MODE="${MODES[$IDX]}"
SEED="${SEEDS[$IDX]}"

if [ "${MODE}" = "index" ]; then
    DV_JSON="outputs/02_calibration/dv_uniformity_calibration/index/calibrated_dv_tolerances.json"
else
    DV_JSON="outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json"
fi

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/wrapper_mode_ablation] array=${IDX} mode=${MODE} seed=${SEED} cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/wrapper_mode_ablation.py \
    --wrapper-mode "${MODE}" \
    --statistic ad \
    --algorithm borg_mm \
    --nfe 200000 \
    --seed "${SEED}" \
    --n-years 20 \
    --ssi 3 \
    --n-islands 1 \
    --checkpoint-freq 10000 \
    --dv-uniformity-json "${DV_JSON}"
