#!/bin/bash
#SBATCH --job-name=04_run_moea_find
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / run_moea_find -- single PARAMETERISED production launcher.
# Replaces the former per-variant slurms (t1_short_block_drb*,
# first_event_ssi3_t10*, primary, ...). Every framing lives as a YAML in
# workflows/04_moea_find_single_site/configs/ and is selected at submit:
#
#   sbatch workflows/04_moea_find_single_site/slurm/run_moea_find.slurm \
#       [CONFIG] [SEED] [SLUG_SUFFIX]
#
#   CONFIG       config basename or path  (default: first_event_ssi3_t10.yaml)
#   SEED         optional --seed override
#   SLUG_SUFFIX  optional --slug-suffix discriminator (e.g. iter1)
#
# 3 nodes x 40 cores = 120-rank MM Borg (exclusive, full rank per policy).
# Calibration prerequisite (must exist for dv_uniform configs):
#   outputs/02_calibration/dv_uniformity_calibration/<mode>/calibrated_dv_tolerances.json

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

CONFIG="${1:-first_event_ssi3_t10.yaml}"
case "${CONFIG}" in
  */*) CONFIG_PATH="${CONFIG}" ;;
  *)   CONFIG_PATH="workflows/04_moea_find_single_site/configs/${CONFIG}" ;;
esac
SEED="${2:-}"
SLUG_SUFFIX="${3:-}"

EXTRA=()
[ -n "${SEED}" ] && EXTRA+=(--seed "${SEED}")
[ -n "${SLUG_SUFFIX}" ] && EXTRA+=(--slug-suffix "${SLUG_SUFFIX}")

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/run_moea_find] config=${CONFIG_PATH} seed=${SEED:-default} \
suffix=${SLUG_SUFFIX:-none} cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python \
    workflows/04_moea_find_single_site/run_moea_find.py \
    --config "${CONFIG_PATH}" "${EXTRA[@]}"
