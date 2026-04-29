# Shared SLURM/bash helpers for MOEA-FIND experiment scripts.
#
# Source this from every .slurm file *before* any work is done:
#     source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
#
# Edit the CLUSTER_* variables below once per HPC and the individual .slurm
# files will inherit them. Nothing in this file is executed at module-load
# time; every command lives inside a function so sourcing is cheap.

# Only enforce strict mode in non-interactive contexts (SLURM jobs, direct execution).
if [[ "${BASH_SOURCE[0]}" == "${0}" ]] || [[ -n "${SLURM_JOB_ID:-}" ]]; then
    set -euo pipefail
fi

# -----------------------------------------------------------------------------
# Cluster-specific configuration — EDIT THESE for your HPC
# -----------------------------------------------------------------------------
CLUSTER_ACCOUNT="${CLUSTER_ACCOUNT:-tja73}"
CLUSTER_PARTITION="${CLUSTER_PARTITION:-normal}"
CLUSTER_PYTHON_MODULE="${CLUSTER_PYTHON_MODULE:-python/3.11.5}"
CLUSTER_MPI_MODULE="${CLUSTER_MPI_MODULE:-openmpi4/4.0.5}"

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# _common.sh lives in workflows/, so PROJECT_ROOT is one dir up.
WORKFLOWS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${WORKFLOWS_DIR}/.." && pwd)"
CLUSTER_VENV="${CLUSTER_VENV:-${PROJECT_ROOT}/venv}"
SLURM_LOG_DIR="${WORKFLOWS_DIR}/slurm/slurm_logs"
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
    # Do NOT module purge — default modules include networking/interconnect
    # components needed for multi-node MPI. Just load what we need.
    module load "${CLUSTER_PYTHON_MODULE}"
}

moea_activate_venv() {
    if [[ -d "${CLUSTER_VENV}" ]]; then
        # shellcheck disable=SC1091
        source "${CLUSTER_VENV}/bin/activate"
    else
        echo "[_common.sh] WARNING: venv not found at ${CLUSTER_VENV}" >&2
    fi
    export PYTHONPATH="${PROJECT_ROOT}:${PROJECT_ROOT}/lib/borg:${PYTHONPATH:-}"
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
    # Borg calls MPI_Init internally so we use mpirun (not srun).
    # --oversubscribe per Water Programming bootcamp pattern.
    # Route OOB and TCP through InfiniBand (ib0) — the Ethernet interface
    # (eno1np0) is firewalled between compute nodes.
    mpirun --oversubscribe -np "${ntasks}" \
        --mca oob_tcp_if_include ib0 \
        --mca btl_tcp_if_include ib0 \
        "$@"
}
