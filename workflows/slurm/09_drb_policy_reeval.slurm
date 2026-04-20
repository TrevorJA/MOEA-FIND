#!/bin/bash
#SBATCH --job-name=moea_09_reeval
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=tja73
#SBATCH --partition=normal
# Manuscript: Fig 9 — Scenario discovery in drought characteristic space
# Experiment: 09 — MOEA-FIND Pareto → Pywr-DRB → metric bank → classifier sweep
# Parallelism: MPI for Stages 2-3 (preprocessing + simulation); serial for 1 & 4
#
# Usage:
#   sbatch workflows/slurm/09_drb_policy_reeval.slurm <path/to/results.json>
#
# Optional environment overrides (set via --export):
#   SUBSET     path to a JSON list of Pareto indices to process (dev ensembles)
#   MANIFEST   path to a satisficing manifest YAML (enables classifier sweep)
#
# Examples:
#   # Full production run, all Pareto members, legacy FFMP classifier only
#   sbatch workflows/slurm/09_drb_policy_reeval.slurm \
#       outputs/exp04_kirsch_single_site/residual_T20_nfe200000_s42_constrained/results.json
#
#   # Dev smoke test: 10 Pareto members + manifest sweep
#   sbatch --export=ALL,SUBSET=outputs/data_cache/dev_10.json,MANIFEST=workflows/experiments/satisficing_manifest.yaml \
#       workflows/slurm/09_drb_policy_reeval.slurm \
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

# Forward optional SUBSET / MANIFEST env vars to each stage invocation.
SUBSET_ARG=()
if [[ -n "${SUBSET:-}" ]]; then
    SUBSET_ARG=(--realization-subset-file "${SUBSET}")
fi

MANIFEST_ARG=()
if [[ -n "${MANIFEST:-}" ]]; then
    MANIFEST_ARG=(--manifest "${MANIFEST}")
fi

echo "[09.slurm] Pareto results: ${PARETO_RESULTS}"
echo "[09.slurm] Subset: ${SUBSET:-<full archive>}"
echo "[09.slurm] Manifest: ${MANIFEST:-<legacy FFMP classifier only>}"
echo "[09.slurm] Processors: ${N_PROCS}"
echo "[09.slurm] Start: $(date)"

# Stage 1 (generate): serial — fits generators, replays DVs
python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage generate \
    --n-years 20 \
    --demand-source constant_max \
    "${SUBSET_ARG[@]}"

# Stages 2-3 (prep + simulate): MPI-parallel
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage prep \
        --n-years 20 \
        --demand-source constant_max \
        --use-mpi \
        "${SUBSET_ARG[@]}"

mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --stage simulate \
        --n-years 20 \
        --demand-source constant_max \
        --batch-size 10 \
        --use-mpi \
        "${SUBSET_ARG[@]}"

# Stage 4 (classify): serial — metric bank, legacy FFMP map, optional sweep
python "${SLURM_SUBMIT_DIR}/workflows/experiments/09_drb_policy_reeval.py" \
    --pareto-results "${PARETO_RESULTS}" \
    --stage classify \
    --n-years 20 \
    --plot \
    "${SUBSET_ARG[@]}" \
    "${MANIFEST_ARG[@]}"

echo "[09.slurm] End: $(date)"
