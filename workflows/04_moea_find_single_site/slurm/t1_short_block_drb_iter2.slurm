#!/bin/bash
#SBATCH --job-name=04_t1_short_block_iter2
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Stage 04 / T=1 short-block iter 2 — anti-ideal re-anchored from iter 1 Pareto.
# D* = 1.5 × Pareto-max for HEADROOM_TIMES_MAX metrics (summer_recession,
# aug_zscore). djf_total_neg / ond_total_neg keep CONSTANT D*=0.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/t1_short_block_iter2] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

REF=outputs/04_moea_find_single_site/run_moea_find/residual_T1_nfe200000_s42_constrained_cmdv_uniform_stad/results.json

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/t1_short_block_drb.yaml \
    --slug-suffix iter2 \
    --anti-ideal-reference "${REF}"
