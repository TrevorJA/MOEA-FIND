#!/bin/bash
#SBATCH --job-name=04_run_moea_find_smoke
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / run_moea_find -- quick wiring SMOKE (replaces the former
# *_smoke variants). Per MM-Borg policy
# (see memory/feedback_mmborg_exclusive_and_full_rank): smokes drop
# NFE, not ranks -- 3 nodes x 40 ranks = 120 cores matching production.
#
#   sbatch workflows/04_moea_find_single_site/slurm/run_moea_find_smoke.slurm \
#       [CONFIG] [NFE]
#
#   CONFIG  config basename/path (default: first_event_ssi3_t10.yaml)
#   NFE     function evals       (default: 10000)

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

CONFIG="${1:-first_event_ssi3_t10.yaml}"
case "${CONFIG}" in
  */*) CONFIG_PATH="${CONFIG}" ;;
  *)   CONFIG_PATH="workflows/04_moea_find_single_site/configs/${CONFIG}" ;;
esac
NFE="${2:-10000}"

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/run_moea_find_smoke] config=${CONFIG_PATH} nfe=${NFE} cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python \
    workflows/04_moea_find_single_site/run_moea_find.py \
    --config "${CONFIG_PATH}" --nfe "${NFE}" --slug-suffix smoke
