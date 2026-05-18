#!/bin/bash
#SBATCH --job-name=04_first_event_ssi3_t10_iter1
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --array=0-3
# Stage 04 / T=10 first-SSI3-drought iter1 production.
# 4 seeds × 200k NFE × 120-rank MM Borg, all running in parallel.
# Anti-ideal anchored on the all-events historical envelope from
# compute_ssi_drought_characteristics on the full 73-year record;
# DV-uniformity AD constraint at cannonsville_T10_ad tolerance 3.169.
# Calibration prerequisites (verified present 2026-05-06):
#   outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json
# Smoke wiring confirmed in job 227084 (Pareto=4455, hyperplane std=0.017).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

SEEDS=(42 7 13 19)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
SEED="${SEEDS[$IDX]}"

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/first_event_ssi3_t10_iter1] array=${IDX} seed=${SEED} cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/first_event_ssi3_t10.yaml \
    --seed "${SEED}" \
    --slug-suffix iter1
