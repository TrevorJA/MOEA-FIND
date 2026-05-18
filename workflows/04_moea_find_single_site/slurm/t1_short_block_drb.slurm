#!/bin/bash
#SBATCH --job-name=04_t1_short_block_drb
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / T=1 short-block event-discovery run (DD-15).
# Metrics: djf_total_neg | summer_recession | aug_zscore | ond_total_neg
# Constraint: DV-uniformity (Anderson-Darling, n_dvs=12).
# Prerequisite: run workflows/02_calibration/slurm/dv_uniformity_T1.slurm first.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/t1_short_block_drb] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/t1_short_block_drb.yaml
