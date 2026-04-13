# Shared SLURM/bash helpers for MOEA-FIND experiment scripts.
#
# Source this from every NN_*.slurm file *before* any work is done:
#     # shellcheck source=./_common.sh
#     source "$(dirname "$0")/_common.sh"
#
# Edit the CLUSTER_* variables below once per HPC and the individual .slurm
# files will inherit them. Nothing in this file is executed at module-load
# time; every command lives inside a function so sourcing is cheap.

set -euo pipefail

# -----------------------------------------------------------------------------
# Cluster-specific configuration — EDIT THESE for your HPC
# -----------------------------------------------------------------------------
CLUSTER_ACCOUNT="${CLUSTER_ACCOUNT:-CHANGE_ME_ACCOUNT}"
CLUSTER_PARTITION="${CLUSTER_PARTITION:-CHANGE_ME_PARTITION}"
CLUSTER_PYTHON_MODULE="${CLUSTER_PYTHON_MODULE:-python/3.11}"
CLUSTER_MPI_MODULE="${CLUSTER_MPI_MODULE:-openmpi/4.1.5}"
CLUSTER_VENV="${CLUSTER_VENV:-${HOME}/.venvs/moea-find}"

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# All .slurm scripts live in scripts/, so PROJECT_ROOT is one dir up.
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"
SLURM_LOG_DIR="${SCRIPTS_DIR}/slurm_logs"
OUTPUT_ROOT="${PROJECT_ROOT}/outputs"
FIGURE_ROOT="${PROJECT_ROOT}/figures"

mkdir -p "${SLURM_LOG_DIR}" "${OUTPUT_ROOT}" "${FIGURE_ROOT}"

# -----------------------------------------------------------------------------
# Environment bring-up
# -----------------------------------------------------------------------------
moea_load_modules() {
    # Make module() available under non-interactive shells.
    if [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    fi
    module purge || true
    module load "${CLUSTER_PYTHON_MODULE}"
    module load "${CLUSTER_MPI_MODULE}"
}

moea_activate_venv() {
    if [[ -d "${CLUSTER_VENV}" ]]; then
        # shellcheck disable=SC1091
        source "${CLUSTER_VENV}/bin/activate"
    else
        echo "[_common.sh] WARNING: venv not found at ${CLUSTER_VENV}" >&2
    fi
    export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
}

moea_setup() {
    moea_load_modules
    moea_activate_venv
    cd "${PROJECT_ROOT}"
    echo "[_common.sh] project=${PROJECT_ROOT}"
    echo "[_common.sh] python=$(command -v python)"
    echo "[_common.sh] host=$(hostname) slurm_job=${SLURM_JOB_ID:-local} task=${SLURM_ARRAY_TASK_ID:-none}"
}

# -----------------------------------------------------------------------------
# MPI launcher selector — prefers srun under SLURM, falls back to mpirun
# -----------------------------------------------------------------------------
moea_mpi_launch() {
    local ntasks="$1"; shift
    if command -v srun >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        srun --ntasks="${ntasks}" --mpi=pmix "$@"
    else
        mpirun -np "${ntasks}" "$@"
    fi
}
