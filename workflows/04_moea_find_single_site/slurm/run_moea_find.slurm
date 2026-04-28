#!/bin/bash
#SBATCH --job-name=moea_04_kirsch_single_site
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# Manuscript: §6.1-6.2 Single-Site Kirsch (Figures 5-6)
# Experiment:  04 — single-site Kirsch MOEA-FIND with MM Borg
# Parallelism: MPI master-worker (MM Borg); 1 master + (N-1) evaluators
#
# Constraint: DV-space uniformity, Anderson-Darling statistic (production default).
# Rationale: AD yields comparable drought-space coverage to the 5-statistic
# hydrologic formulation with a tighter, more interpretable constraint geometry.
# Calibration: outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json
# (produced by diag_dv_uniformity_calibration.slurm). Must exist before submitting.
#
# Defaults (NFE, n_years, seed, mode, ssi timescale, islands, checkpoint
# frequency) come from src/experiment_config.py. Override per-run via env vars:
#
#   sbatch --export=ALL,BORG_NFE=200000 workflows/04_moea_find_single_site/slurm/run_moea_find.slurm
#   sbatch --export=ALL,BORG_NFE=200000,BORG_SEED=101 workflows/04_moea_find_single_site/slurm/run_moea_find.slurm
#
# Smoke test (shrink allocation):
#   sbatch --export=ALL,BORG_NFE=2000 --nodes=1 --ntasks-per-node=8 --time=00:15:00 \
#          workflows/04_moea_find_single_site/slurm/run_moea_find.slurm
#
# To reproduce the SI hydrologic-constraint arm instead:
#   sbatch --export=ALL,BORG_NFE=200000,CONSTRAINT_MODE=hydrologic \
#          workflows/04_moea_find_single_site/slurm/run_moea_find.slurm

# SLURM copies scripts to a spool dir; use SLURM_SUBMIT_DIR instead.
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

# Build the Python arg list, passing only CLI overrides that were set.
# Omitted flags fall back to defaults in src/experiment_config.py.
# Use ${VAR:-} (colon-dash) because _common.sh enables `set -u`: unbound
# variables would otherwise error instead of being treated as empty.
# Production constraint defaults — override via --export if needed.
: "${CONSTRAINT_MODE:=dv_uniform}"
: "${DV_STATISTIC:=ad}"
: "${DV_UNIFORMITY_JSON:=outputs/diag_dv_uniformity_calibration/calibrated_dv_tolerances.json}"

PY_ARGS=(--algorithm borg_mm --plot)
PY_ARGS+=(--constraint-mode "${CONSTRAINT_MODE}")
PY_ARGS+=(--statistic       "${DV_STATISTIC}")
PY_ARGS+=(--dv-uniformity-json "${DV_UNIFORMITY_JSON}")
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
moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    "${PY_ARGS[@]}"

# Package outputs to SCRATCH for safekeeping
RESULTS_DIR="${PROJECT_ROOT}/outputs/exp04_kirsch_single_site"
SCRATCH_DEST="${SCRATCH:-/tmp}/moea_find/exp04_${SLURM_JOB_ID}"
mkdir -p "${SCRATCH_DEST}"
cp -r "${RESULTS_DIR}" "${SCRATCH_DEST}/" 2>/dev/null || true
echo "[slurm] packaged results to ${SCRATCH_DEST}"
