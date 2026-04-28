#!/bin/bash
#SBATCH --job-name=moea_09_reeval
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=30
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal
# MPI transport note (Hopper, OpenMPI 4.0.5):
# The original flags `--mca oob/btl_tcp_if_include ib0` did NOT prevent
# OpenMPI from attempting to load the `openib` BTL (RDMA), and that
# attempt hangs indefinitely on this fabric — see diagnostic jobs
# 212674/212678 (2026-04-22) which isolated this to `MPI_ERR_OTHER` at
# the first collective whenever the `openib` BTL is in play.  The fix
# (applied below on each mpirun) is an explicit PML/BTL whitelist:
#     --mca pml ob1 --mca btl self,vader,tcp
# which uses only local (self), shared-memory (vader), and TCP — never
# openib — and passes a 90-rank multi-node Barrier in ~3 seconds.
# Experiment: 09 — MOEA-FIND Pareto → Pywr-DRB → metric bank + legacy
# satisficing table. Scenario-discovery figures (GBT / logreg multi-panel
# maps, AUC summaries, legacy single-panel map) are emitted by
# workflows/07_scenario_discovery/scenario_discovery_plots.py, which reads
# this script's persisted outputs.
# Parallelism: MPI for Stages 2-3 (preprocessing + simulation); serial for 1 & 4
#
# Usage:
#   sbatch workflows/06_pywrdrb_reeval/slurm/policy_reeval.slurm <path/to/results.json>
#
# Optional environment overrides (set via --export):
#   SUBSET          path to a JSON list of Pareto indices to process (dev ensembles)
#   SKIP_GENERATE=1 skip Stage 1 (resume from Stage 2 when Stage 1 output exists)
#   SKIP_PREP=1     skip Stage 2 (resume from Stage 3 when predicted_inflows_mgd.hdf5 exists)
#
# Examples:
#   # Full production run, all Pareto members
#   sbatch workflows/06_pywrdrb_reeval/slurm/policy_reeval.slurm \
#       outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json
#
#   # Dev smoke test: 10 Pareto members
#   sbatch --export=ALL,SUBSET=outputs/data_cache/dev_10.json \
#       workflows/06_pywrdrb_reeval/slurm/policy_reeval.slurm \
#       outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Unbuffered stdout so progress prints from long-running loops
# (e.g. per-realization replay in Stage 1) land in the SLURM log
# in real time rather than only on process exit.
export PYTHONUNBUFFERED=1

# Pin every BLAS/OpenMP backend to a single thread. Stage 1 uses
# joblib to parallelise Nowak disaggregation across 16 worker
# processes; without this, each worker's numpy/sklearn internals
# try to use all 16 cores simultaneously, so 16 workers × 16 threads
# thrash the CPU scheduler and we see ~1x speedup instead of ~16x.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

PARETO_RESULTS="${1:?Usage: sbatch 09_drb_policy_reeval.slurm <path/to/results.json>}"
N_PROCS="${SLURM_NTASKS:-16}"

# Forward optional SUBSET env var to each stage invocation.
SUBSET_ARG=()
if [[ -n "${SUBSET:-}" ]]; then
    SUBSET_ARG=(--realization-subset-file "${SUBSET}")
fi

echo "[09.slurm] Pareto results: ${PARETO_RESULTS}"
echo "[09.slurm] Subset: ${SUBSET:-<full archive>}"
echo "[09.slurm] Processors: ${N_PROCS}"
echo "[09.slurm] Start: $(date)"

# Stage 1 (generate): serial — fits generators, replays DVs
if [[ "${SKIP_GENERATE:-0}" == "1" ]]; then
    echo "[09.slurm] Skipping Stage 1 (SKIP_GENERATE=1)"
else
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage generate \
        --n-years 20 \
        --demand-source constant_max \
        "${SUBSET_ARG[@]}"
fi

# Stages 2-3 (prep + simulate): MPI-parallel
if [[ "${SKIP_PREP:-0}" == "1" ]]; then
    echo "[09.slurm] Skipping Stage 2 (SKIP_PREP=1)"
else
    mpirun --oversubscribe -np "${N_PROCS}" \
        --mca pml ob1 \
        --mca btl self,vader,tcp \
        python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
            --pareto-results "${PARETO_RESULTS}" \
            --stage prep \
            --n-years 20 \
            --demand-source constant_max \
            --use-mpi \
            "${SUBSET_ARG[@]}"
fi

mpirun --oversubscribe -np "${N_PROCS}" \
    --mca pml ob1 \
    --mca btl self,vader,tcp \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage simulate \
        --n-years 20 \
        --demand-source constant_max \
        --batch-size 5 \
        --use-mpi \
        "${SUBSET_ARG[@]}"

# Stage 4 (classify): serial — metric bank + legacy FFMP classification
# table.  Downstream: run 17_drb_scenario_discovery_plots.slurm for
# figures (single-panel, multi-panel GBT/logreg, AUC summaries).
python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage classify \
    --n-years 20 \
    "${SUBSET_ARG[@]}"

echo "[09.slurm] End: $(date)"
