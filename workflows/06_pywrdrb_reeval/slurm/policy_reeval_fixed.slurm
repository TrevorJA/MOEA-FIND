#!/bin/bash
#SBATCH --job-name=06_policy_reeval_fixed
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal

# Stage 06 / policy_reeval -- MOEA-FIND Pareto -> Pywr-DRB -> metric bank.
# Manuscript: contributes data to Fig 9 (scenario discovery, stage 07).
#
# MPI transport note (Hopper, OpenMPI 4.0.5):
# The original flags `--mca oob/btl_tcp_if_include ib0` did NOT prevent
# OpenMPI from attempting to load the `openib` BTL (RDMA), and that
# attempt hangs indefinitely on this fabric -- see diagnostic jobs
# 212674/212678 (2026-04-22) which isolated this to MPI_ERR_OTHER at
# the first collective whenever the `openib` BTL is in play. The fix
# (applied below on each mpirun) is an explicit PML/BTL whitelist:
#     --mca pml ob1 --mca btl self,vader,tcp
# which uses only local (self), shared-memory (vader), and TCP -- never
# openib -- and passes a 90-rank multi-node Barrier in ~3 seconds.
#
# Parallelism: 10 nodes x 30 ntasks-per-node = 300 ranks for stage 3
# (simulation). Stages 1 + 4 run serially on rank 0. Full HPC budget.
#
# History of configurations tried on Hopper (April 2026):
#   - 3x30 = 90 ranks (no --exclusive): OOM-killed at 3.5 h (cgroup memory
#     cap + memory growth in pywrdrb sims).
#   - 5x40 = 200 ranks --exclusive: OOM in 4 min (40 ranks/node x ~5 GB
#     per rank exceeded Hopper's ~192 GB/node).
#   - 8x20 = 160 ranks --exclusive: 100 ranks came up healthy and wrote sims;
#     ranks 100-159 died silently at MPI startup; MPI Barrier in
#     pywrdrb_bridge.py:758 hung waiting on the dead 60.
#
# Two more failed attempts after the 5x40/8x20 attempts above:
#   - 10x30 = 300 ranks --exclusive, batch_size=10: OOM at 7 min (rank 29)
#   - 10x30 = 300 ranks --exclusive, batch_size=5:  OOM at 18 min (rank 27)
# 30 ranks/node is over the per-rank memory cap on Hopper compute nodes
# even with --exclusive. The 18-min attempt also showed an InfiniBand
# connect error on c0009 ("connect() to 192.168.12.10:1031 failed").
#
# Current (10x20 = 200 ranks --exclusive): each node hosts only 20 ranks
# with ~9.5 GB per-rank memory headroom -- the same density that ran
# healthy at 5x20=100. 3308 sims / 200 ranks ~= 17 sims per rank = 4 batches
# of 5, ~70 min wall for simulate. batch_size=5 is the production cap.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Bake the Pareto-archive path. To re-evaluate a different stage 04 run,
# edit this line; do NOT pass arguments via sbatch.
PARETO_RESULTS="outputs/04_moea_find_single_site/run_moea_find/residual_T10_nfe200000_s42_constrained_cmdv_uniform_sfxfixed_stad/results.json"

# Unbuffered stdout so progress prints from long-running loops
# (e.g. per-realization replay in Stage 1) land in the SLURM log
# in real time rather than only on process exit.
export PYTHONUNBUFFERED=1

# Pin every BLAS/OpenMP backend to a single thread. Stage 1 uses
# joblib to parallelise Nowak disaggregation across 16 worker
# processes; without this, each worker's numpy/sklearn internals
# try to use all 16 cores simultaneously, so 16 workers x 16 threads
# thrash the CPU scheduler and we see ~1x speedup instead of ~16x.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

N_PROCS="${SLURM_NTASKS:-200}"

echo "[06.slurm] Pareto results: ${PARETO_RESULTS}"
echo "[06.slurm] Processors: ${N_PROCS}"
echo "[06.slurm] Start: $(date)"

# --- Stage 1 (generate) + Stage 2 (prep): DONE SEPARATELY. ---
# generate: SKIPPED (outputs persist for this NEW slug: gage_flow_mgd.hdf5,
#   catchment_inflow_mgd.hdf5; verified on disk).
# prep: moved to its own SINGLE-NODE job
#   workflows/06_pywrdrb_reeval/slurm/policy_reeval_prep_fixed.slurm
#   because Hopper inter-node IB/TCP is unreliable and even rank-placement
#   isolation left ORTE/OOB spanning all 10 nodes (empty presim HDF5, all
#   reservoirs missing). This script REQUIRES prep's outputs already on disk:
#   pywrdrb_inputs/{predicted_inflows_mgd.hdf5, presimulated_releases_mgd.hdf5}.
# This job runs ONLY Stage 3 (simulate, memory-tuned 10-node/200-rank) +
# Stage 4 (classify, serial).

# Stage 3 (simulate): MPI-parallel. Per-rank memory pressure is now
# managed by --exclusive + 5-node allocation (full node memory available
# per rank), keeping batch_size=5 as the production default.
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca pml ob1 \
    --mca btl self,vader,tcp \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage simulate \
        --n-years 10 \
        --demand-source constant_max \
        --batch-size 5 \
        --use-mpi

# Stage 4 (classify): serial -- metric bank + legacy FFMP table.
python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage classify \
    --n-years 10

echo "[06.slurm] End: $(date)"
