#!/bin/bash
#SBATCH --job-name=moea_15_wrapper_mode_ablation
#SBATCH --output=workflows/slurm/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --array=0-5
# Manuscript Supporting Information: Kirsch wrapper-mode ablation (script 15).
# Array index maps to (mode, seed):
#   0..2 -> index x {42, 101, 2026}
#   3..5 -> residual x {42, 101, 2026}
# Each task: independent MM Borg MPI run, identical NFE / DV length / constraint.
# Both modes use DV-uniformity AD constraint.
#
# Precondition (calibration must exist before submitting):
#   outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json
#         (sbatch workflows/slurm/diag_dv_uniformity_calibration.slurm)
#
# Usage:
#   sbatch --export=ALL,BORG_NFE=200000 workflows/slurm/15_wrapper_mode_ablation.slurm
# Tight budget (smaller NFE):
#   sbatch --export=ALL,BORG_NFE=50000 workflows/slurm/15_wrapper_mode_ablation.slurm
# Single mode/seed for a smoke test:
#   sbatch --export=ALL,BORG_NFE=2000 --array=0 \
#          --time=00:30:00 --nodes=1 --ntasks-per-node=8 \
#          workflows/slurm/15_wrapper_mode_ablation.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

MODES=(index index index residual residual residual)
SEEDS=(42 101 2026 42 101 2026)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
MODE="${MODES[$IDX]}"
SEED="${SEEDS[$IDX]}"

: "${BORG_NFE:=200000}"
: "${N_YEARS:=20}"
: "${SSI_TIMESCALE:=3}"
: "${N_ISLANDS:=1}"
: "${CHECKPOINT_FREQ:=10000}"

PY_ARGS=(
    --wrapper-mode "${MODE}"
    --statistic ad
    --algorithm borg_mm
    --nfe "${BORG_NFE}"
    --seed "${SEED}"
    --n-years "${N_YEARS}"
    --ssi "${SSI_TIMESCALE}"
    --n-islands "${N_ISLANDS}"
    --checkpoint-freq "${CHECKPOINT_FREQ}"
)

# The default DV-uniformity calibration JSON has n_dvs=240 (residual mode,
# n_years_out*12). Index mode needs n_dvs=252 ((n_years_out+1)*12) because
# SynHydro's generate_from_indices requires an extra year buffer for Dec-Jan
# handling. Point index-mode runs at a separate, mode-specific calibration.
if [ "${MODE}" = "index" ]; then
    PY_ARGS+=(--dv-uniformity-json outputs/diag_dv_uniformity_calibration_index/calibrated_dv_tolerances.json)
fi

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[slurm] array=${IDX} mode=${MODE} seed=${SEED} nfe=${BORG_NFE} "
echo "[slurm] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"
echo "[slurm] python args: ${PY_ARGS[*]}"
moea_mpi_launch "${N_PROCESSORS}" python workflows/experiments/15_wrapper_mode_ablation.py \
    "${PY_ARGS[@]}"

# Package per-mode results to scratch so they survive run_dir rotation.
RESULTS_DIR="${PROJECT_ROOT}/outputs/exp15_wrapper_mode_ablation/${MODE}"
SCRATCH_DEST="${SCRATCH:-/tmp}/moea_find/exp15_${SLURM_ARRAY_JOB_ID:-local}_${MODE}_s${SEED}"
mkdir -p "${SCRATCH_DEST}"
cp -r "${RESULTS_DIR}" "${SCRATCH_DEST}/" 2>/dev/null || true
echo "[slurm] packaged ${RESULTS_DIR} -> ${SCRATCH_DEST}"
