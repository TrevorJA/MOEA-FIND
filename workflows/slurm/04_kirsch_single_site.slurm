#!/bin/bash
#SBATCH --job-name=moea_04_kirsch_single_site
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Manuscript: §6.1-6.2 Single-Site Kirsch (Figures 5-6)
# Experiment:  04 — single-site Kirsch MOEA-FIND with MM Borg
# Parallelism: MPI master-worker (MM Borg); 1 master + (N-1) evaluators
#
# IMPORTANT: Run diag_constraint_calibration.slurm FIRST to produce
# calibrated_tolerances.json. Without it the run proceeds unconstrained.
#
# Override variables in the sbatch command line:
#   sbatch --export=ALL,BORG_NFE=200000,BORG_MODE=index scripts/04_kirsch_single_site.slurm

# SLURM copies scripts to a spool dir; use SLURM_SUBMIT_DIR instead.
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

: "${BORG_NFE:=500000}"
: "${BORG_MODE:=residual}"
: "${BORG_SEED:=42}"
: "${N_YEARS:=20}"
: "${SSI_TIMESCALE:=3}"
: "${N_ISLANDS:=1}"
: "${CHECKPOINT_FREQ:=50000}"

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

# Path to calibration JSON produced by diag_constraint_calibration.slurm.
CONSTRAINTS_JSON="${PROJECT_ROOT}/outputs/diag_constraint_calibration/calibrated_tolerances.json"

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[slurm] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"
moea_mpi_launch "${N_PROCESSORS}" python workflows/experiments/04_kirsch_single_site.py \
    --nfe "${BORG_NFE}" \
    --mode "${BORG_MODE}" \
    --n-years "${N_YEARS}" \
    --ssi "${SSI_TIMESCALE}" \
    --seed "${BORG_SEED}" \
    --algorithm borg_mm \
    --constraints-json "${CONSTRAINTS_JSON}" \
    --n-islands "${N_ISLANDS}" \
    --checkpoint-freq "${CHECKPOINT_FREQ}" \
    --plot

# Package outputs to SCRATCH for safekeeping
RESULTS_DIR="${PROJECT_ROOT}/outputs/exp04_kirsch_single_site"
SCRATCH_DEST="${SCRATCH:-/tmp}/moea_find/exp04_${SLURM_JOB_ID}"
mkdir -p "${SCRATCH_DEST}"
cp -r "${RESULTS_DIR}" "${SCRATCH_DEST}/" 2>/dev/null || true
echo "[slurm] packaged results to ${SCRATCH_DEST}"
