#!/bin/bash
#SBATCH --job-name=04_t1_short_block_drb_v2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / T=1 short-block run — DD-15c reformulation (bounded metrics).
# Metrics: oct_logmean_g | dec_logmean_g | jun_logmean_g | recession_amjjas_g
# All bounded in [0, 1) via Mapping G; anti-ideal D* = 1.0 (CONSTANT).
# Constraint: DV-uniformity (Anderson-Darling, n_dvs=12).
# Compares against the v1 baseline at:
#   outputs/04_moea_find_single_site/run_moea_find/residual_T1_nfe200000_s42_constrained_cmdv_uniform_stad/

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/t1_short_block_drb_v2] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/t1_short_block_drb_v2.yaml \
    --slug-suffix v2
