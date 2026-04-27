#!/bin/bash
#SBATCH --job-name=diag_mpi_single_node
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --exclusive
#SBATCH --time=00:02:00
#SBATCH --account=tja73
#SBATCH --partition=normal
#SBATCH --exclude=c[0001-0003]
# Single-node MPI sanity: if multi-node Barrier is broken but single-node works,
# the issue is the inter-node transport (IB/TCP OOB), not the MPI stack itself.
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup
export PYTHONUNBUFFERED=1

echo "[single-node] ranks=${SLURM_NTASKS} nodelist=${SLURM_JOB_NODELIST}"
mpirun --oversubscribe -np "${SLURM_NTASKS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python -c "
from mpi4py import MPI
import socket
comm = MPI.COMM_WORLD
comm.Barrier()
if comm.Get_rank() == 0:
    print(f'single-node barrier ok size={comm.Get_size()}', flush=True)
print(f'rank {comm.Get_rank()} host {socket.gethostname()}', flush=True)
"
echo "[single-node] rc=$?"
