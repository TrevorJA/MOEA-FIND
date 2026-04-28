#!/bin/bash
#SBATCH --job-name=diag_mpi_variant
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=30
#SBATCH --exclusive
#SBATCH --time=00:03:00
#SBATCH --account=tja73
#SBATCH --partition=normal
#SBATCH --exclude=c[0001-0003]
# Usage:
#   sbatch --export=ALL,VARIANT=B diag_mpi_variant.slurm   # mpirun no flags
#   sbatch --export=ALL,VARIANT=C diag_mpi_variant.slurm   # srun --mpi=pmix
#   sbatch --export=ALL,VARIANT=A diag_mpi_variant.slurm   # prod flags (expected to hang)
#   sbatch --export=ALL,VARIANT=E diag_mpi_variant.slurm   # mpirun -mca pml ucx (UCX transport)

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1
N_PROCS="${SLURM_NTASKS:-90}"
V="${VARIANT:?set VARIANT=A|B|C|E via --export}"
PY='
import sys, socket
from mpi4py import MPI
comm = MPI.COMM_WORLD
r, sz = comm.Get_rank(), comm.Get_size()
if r == 0:
    print(f"pre-barrier: size={sz}", flush=True)
comm.Barrier()
if r == 0:
    print(f"barrier1 ok", flush=True)
data = comm.bcast("hello" if r == 0 else None, root=0)
comm.Barrier()
if r == 0:
    print(f"bcast+barrier2 ok data={data!r}", flush=True)
print(f"rank {r}/{sz} on {socket.gethostname()}", flush=True)
'

echo "[diag ${V}] ranks=${N_PROCS} nodelist=${SLURM_JOB_NODELIST} start=$(date +%H:%M:%S.%N)"

case "${V}" in
    A)
        mpirun --oversubscribe -np "${N_PROCS}" \
            --mca oob_tcp_if_include ib0 \
            --mca btl_tcp_if_include ib0 \
            python -c "${PY}"
        ;;
    B)
        mpirun --oversubscribe -np "${N_PROCS}" python -c "${PY}"
        ;;
    C)
        srun --mpi=pmix -n "${N_PROCS}" python -c "${PY}"
        ;;
    E)
        mpirun --oversubscribe -np "${N_PROCS}" \
            --mca pml ucx --mca btl ^openib,tcp \
            python -c "${PY}"
        ;;
    G)
        # srun with PMI-2 (the actually-available PMI on this SLURM build)
        srun --mpi=pmi2 -n "${N_PROCS}" python -c "${PY}"
        ;;
    H)
        # mpirun but explicit BTLs, no OOB-over-IB coercion
        mpirun --oversubscribe -np "${N_PROCS}" \
            --mca pml ob1 --mca btl self,vader,tcp \
            python -c "${PY}"
        ;;
    *) echo "unknown VARIANT=${V}"; exit 2 ;;
esac
rc=$?
echo "[diag ${V}] rc=${rc} end=$(date +%H:%M:%S.%N)"
