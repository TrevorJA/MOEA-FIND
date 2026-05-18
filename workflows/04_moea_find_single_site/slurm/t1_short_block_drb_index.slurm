#!/bin/bash
#SBATCH --job-name=04_t1_short_block_index
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / T=1 short-block A/B with INDEX mode.
# Index mode bounds flows to historical residuals by construction (DVs map
# to year indices, not residual quantiles). Tests whether the implausibly
# extreme single-month flows seen in the residual baseline (max 33,020 cfs
# vs historical 7,700) are a residual-tail-sampling artifact.
# Prerequisite: outputs/02_calibration/dv_uniformity_calibration/index/calibrated_dv_tolerances.json
# must contain cannonsville_T1_ad (n_dvs=24).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/t1_short_block_index] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/t1_short_block_drb.yaml \
    --mode index \
    --dv-uniformity-json outputs/02_calibration/dv_uniformity_calibration/index/calibrated_dv_tolerances.json
