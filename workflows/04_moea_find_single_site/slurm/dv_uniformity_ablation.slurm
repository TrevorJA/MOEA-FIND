#!/bin/bash
#SBATCH --job-name=04_dv_uniformity_ablation
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --array=0-5%3
# Stage 04 / dv_uniformity_ablation -- compute only.
# SI-G. 3 arms x 2 seeds = 6 tasks; throttle %3 caps at 360 cores
# (well under the 400-core 10-node x 40-core budget).
# Array index -> (arm, statistic, seed):
#   0: hydrologic       (l2_star ignored)  s=42
#   1: hydrologic       (l2_star ignored)  s=101
#   2: dv_uniform       l2_star            s=42
#   3: dv_uniform       l2_star            s=101
#   4: dv_uniform       ad                 s=42
#   5: dv_uniform       ad                 s=101
# Calibration prerequisites (must exist):
#   outputs/02_calibration/constraint_calibration/calibrated_tolerances.json
#   outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

ARMS=(hydrologic hydrologic dv_uniform dv_uniform dv_uniform dv_uniform)
SEEDS=(42 101 42 101 42 101)
STATISTICS=(l2_star l2_star l2_star l2_star ad ad)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
ARM="${ARMS[$IDX]}"
SEED="${SEEDS[$IDX]}"
STATISTIC="${STATISTICS[$IDX]}"

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/dv_uniformity_ablation] array=${IDX} arm=${ARM} seed=${SEED} stat=${STATISTIC} cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/dv_uniformity_ablation.py \
    --arm "${ARM}" \
    --statistic "${STATISTIC}" \
    --algorithm borg_mm \
    --nfe 200000 \
    --seed "${SEED}" \
    --mode residual \
    --n-years 20 \
    --ssi 3 \
    --n-islands 1 \
    --checkpoint-freq 10000 \
    --constraints-json outputs/02_calibration/constraint_calibration/calibrated_tolerances.json \
    --dv-uniformity-json outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json
