#!/bin/bash
#SBATCH --job-name=moea_13_dv_uniformity_ablation
#SBATCH --output=workflows/slurm/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --array=0-8
# Manuscript Supporting Information: DV-uniformity ablation (script 13).
# Array index maps to (arm, seed):
#   0..2 -> hydrologic x {42, 101, 2026}
#   3..5 -> dv_uniform (l2_star) x {42, 101, 2026}
#   6..8 -> dv_uniform (ad) x {42, 101, 2026}
# Each task: independent MM Borg MPI run, identical NFE / DV length / mode.
#
# Preconditions (both calibrations must exist before submitting):
#   1. outputs/diag_constraint_calibration/calibrated_tolerances.json
#         (sbatch workflows/slurm/diag_constraint_calibration.slurm)
#   2. outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json
#         (sbatch workflows/slurm/diag_dv_uniformity_calibration.slurm)
#
# Usage:
#   sbatch --export=ALL,BORG_NFE=500000 workflows/slurm/13_dv_uniformity_ablation.slurm
# Tight budget (smaller NFE):
#   sbatch --export=ALL,BORG_NFE=250000 workflows/slurm/13_dv_uniformity_ablation.slurm
# Single seed for a smoke test:
#   sbatch --export=ALL,BORG_NFE=2000 --array=0 \
#          --time=00:30:00 --nodes=1 --ntasks-per-node=8 \
#          workflows/slurm/13_dv_uniformity_ablation.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

ARMS=(hydrologic hydrologic hydrologic dv_uniform dv_uniform dv_uniform dv_uniform dv_uniform dv_uniform)
SEEDS=(42 101 2026 42 101 2026 42 101 2026)
STATISTICS=(l2_star l2_star l2_star l2_star l2_star l2_star ad ad ad)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
ARM="${ARMS[$IDX]}"
SEED="${SEEDS[$IDX]}"
STATISTIC="${STATISTICS[$IDX]}"

: "${BORG_NFE:=200000}"
: "${BORG_MODE:=residual}"
: "${N_YEARS:=20}"
: "${SSI_TIMESCALE:=3}"
: "${N_ISLANDS:=1}"
: "${CHECKPOINT_FREQ:=10000}"

PY_ARGS=(
    --arm "${ARM}"
    --statistic "${STATISTIC}"
    --algorithm borg_mm
    --nfe "${BORG_NFE}"
    --seed "${SEED}"
    --mode "${BORG_MODE}"
    --n-years "${N_YEARS}"
    --ssi "${SSI_TIMESCALE}"
    --n-islands "${N_ISLANDS}"
    --checkpoint-freq "${CHECKPOINT_FREQ}"
)

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[slurm] array=${IDX} arm=${ARM} seed=${SEED} nfe=${BORG_NFE} "
echo "[slurm] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"
echo "[slurm] python args: ${PY_ARGS[*]}"
moea_mpi_launch "${N_PROCESSORS}" python workflows/experiments/13_dv_uniformity_ablation.py \
    "${PY_ARGS[@]}"

# Package per-arm results to scratch so they survive run_dir rotation.
RESULTS_DIR="${PROJECT_ROOT}/outputs/exp13_dv_uniformity_ablation/${ARM}"
SCRATCH_DEST="${SCRATCH:-/tmp}/moea_find/exp13_${SLURM_ARRAY_JOB_ID:-local}_${ARM}_s${SEED}"
mkdir -p "${SCRATCH_DEST}"
cp -r "${RESULTS_DIR}" "${SCRATCH_DEST}/" 2>/dev/null || true
echo "[slurm] packaged ${RESULTS_DIR} -> ${SCRATCH_DEST}"
