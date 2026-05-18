#!/bin/bash
#SBATCH --job-name=09_run_mv_sa_coarse
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=19
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --exclusive
#SBATCH --account=tja73
#SBATCH --partition=normal
# Single-node, exclusive: 19 ranks each pinned to 4 cores (with thread
# caps below preventing NumPy fanout) on a node with no other tenants.
# 19 x 4 = 76 cores reserved out of ~80 available; the 4-core slot per
# rank gives memory-bandwidth headroom for N=3308 where SALib delta is
# O(N^2) per Plischke resample. Multi-node MPI on Hopper hangs at
# startup (job 217414); single-node exclusive sidesteps that.
# Stage 09 / run_mv_sa -- magnitude-varying sensitivity analysis (MPI).
# Manuscript: §3.4 / §4.x (paired with §3.3 scenario discovery).
# Pair with: workflows/09_magnitude_varying_sa/slurm/plots/run_mv_sa.slurm
#
# 19 MPI ranks, one per percentile in the default {0.05, 0.10, ..., 0.95}
# grid. Each rank computes its assigned percentile slices independently;
# rank 0 gathers and writes parquet. Bank/chars paths are baked: change
# them by editing this file, not via --export.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

# Thread caps: every MPI rank holds one core (--cpus-per-task=1) and
# we must prevent NumPy/SciPy/SALib from spawning OpenBLAS/MKL/OpenMP
# threads underneath, which would oversubscribe the node by ~19x and
# slow per-rank wall by an order of magnitude. Single-threaded ranks
# are the right choice when the parallelism is across percentiles.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

moea_mpi_launch 19 python workflows/09_magnitude_varying_sa/run_mv_sa.py \
    --bank  outputs/06_pywrdrb_reeval/policy_reeval/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxcoarse_stad/results/metric_bank.csv \
    --chars outputs/06_pywrdrb_reeval/policy_reeval/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxcoarse_stad/results/results_sa.json \
    --magnitude-axis nyc_min_storage_frac \
    --config workflows/09_magnitude_varying_sa/configs/delta_within_trace.yaml
