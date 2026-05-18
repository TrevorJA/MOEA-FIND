#!/bin/bash
#SBATCH --job-name=04_run_moea_find
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / run_moea_find -- compute only.
# Default invocation: T=10 first-SSI3-drought, 5-objective preset
# first_event_ssi3_t10, single seed. For the production seed sweep use
# first_event_ssi3_t10_iter1.slurm. Legacy 3-objective T=20 primary preset
# remains available via --config workflows/04_moea_find_single_site/configs/primary.yaml.
# 3 nodes x 40 cores = 120 cores.
# Calibration prerequisite (must exist):
#   outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/run_moea_find] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/first_event_ssi3_t10.yaml
