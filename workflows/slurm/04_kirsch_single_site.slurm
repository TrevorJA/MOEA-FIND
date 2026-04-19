#!/bin/bash
#SBATCH --job-name=moea_04_kirsch_single_site
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Manuscript: §6.1-6.2 Single-Site Kirsch (Figures 5-6)
# Experiment:  04 — single-site Kirsch MOEA-FIND with MM Borg
# Parallelism: MPI master-worker (MM Borg); 1 master + (N-1) evaluators
#
# Defaults (NFE, n_years, seed, mode, ssi timescale, islands, checkpoint
# frequency, constraints JSON) come from src/experiment_config.py.
# Override per-run from sbatch via env vars:
#
#   sbatch --export=ALL,BORG_NFE=500000 workflows/slurm/04_kirsch_single_site.slurm
#   sbatch --export=ALL,BORG_NFE=100,BORG_MODE=index workflows/slurm/04_kirsch_single_site.slurm
#
# For a tiny smoke test, also shrink the allocation:
#   sbatch --export=ALL,BORG_NFE=100 --nodes=1 --ntasks-per-node=8 --time=00:15:00 \
#          workflows/slurm/04_kirsch_single_site.slurm
#
# IMPORTANT: Run diag_constraint_calibration.slurm FIRST to produce
# calibrated_tolerances.json. Without it the run proceeds unconstrained.

# SLURM copies scripts to a spool dir; use SLURM_SUBMIT_DIR instead.
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

# Build the Python arg list, passing only CLI overrides that were set.
# Omitted flags fall back to defaults in src/experiment_config.py.
# Use ${VAR:-} (colon-dash) because _common.sh enables `set -u`: unbound
# variables would otherwise error instead of being treated as empty.
PY_ARGS=(--algorithm borg_mm --plot)
[ -n "${BORG_NFE:-}"         ] && PY_ARGS+=(--nfe             "${BORG_NFE}")
[ -n "${BORG_MODE:-}"        ] && PY_ARGS+=(--mode            "${BORG_MODE}")
[ -n "${BORG_SEED:-}"        ] && PY_ARGS+=(--seed            "${BORG_SEED}")
[ -n "${N_YEARS:-}"          ] && PY_ARGS+=(--n-years         "${N_YEARS}")
[ -n "${SSI_TIMESCALE:-}"    ] && PY_ARGS+=(--ssi             "${SSI_TIMESCALE}")
[ -n "${N_ISLANDS:-}"        ] && PY_ARGS+=(--n-islands       "${N_ISLANDS}")
[ -n "${CHECKPOINT_FREQ:-}"  ] && PY_ARGS+=(--checkpoint-freq "${CHECKPOINT_FREQ}")
[ -n "${CONSTRAINTS_JSON:-}" ] && PY_ARGS+=(--constraints-json "${CONSTRAINTS_JSON}")

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[slurm] nodes=${SLURM_NNODES} tasks_per_node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"
echo "[slurm] python args: ${PY_ARGS[*]}"
moea_mpi_launch "${N_PROCESSORS}" python workflows/experiments/04_kirsch_single_site.py \
    "${PY_ARGS[@]}"

# Package outputs to SCRATCH for safekeeping
RESULTS_DIR="${PROJECT_ROOT}/outputs/exp04_kirsch_single_site"
SCRATCH_DEST="${SCRATCH:-/tmp}/moea_find/exp04_${SLURM_JOB_ID}"
mkdir -p "${SCRATCH_DEST}"
cp -r "${RESULTS_DIR}" "${SCRATCH_DEST}/" 2>/dev/null || true
echo "[slurm] packaged results to ${SCRATCH_DEST}"
