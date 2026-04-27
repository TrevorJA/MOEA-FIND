#!/bin/bash
#SBATCH --job-name=diag_mpi_bootstrap
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=30
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --account=tja73
#SBATCH --partition=normal
#SBATCH --exclude=c[0001-0003]
# Minimal reproducer for the Stage-3-cold-start Barrier failure seen in
# 212524/212576/212618/212650. Strips everything down to a bare MPI
# Barrier launched with the same mpirun flags the 09 slurm uses.
# Goal: determine whether the MPI stack is the culprit (this job fails)
# or something in our pywrdrb init path is (this job succeeds, but the
# full 09 Stage 3 fails).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export PYTHONUNBUFFERED=1

N_PROCS="${SLURM_NTASKS:-90}"
echo "[diag] ranks=${N_PROCS} nodes=${SLURM_JOB_NUM_NODES} nodelist=${SLURM_JOB_NODELIST}"
echo "[diag] module list:"
module list 2>&1 | sed 's/^/[diag]   /'
echo "[diag] start: $(date +%H:%M:%S.%N)"

# ----- Variant A: mpirun + same MCA flags as 09 (expected to reproduce)
echo
echo "=============================================================="
echo "[diag A] mpirun --oversubscribe -np N --mca oob/btl ib0   (09 prod flags)"
echo "=============================================================="
t0=$(date +%s)
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python -c "
from mpi4py import MPI
import socket, sys
comm = MPI.COMM_WORLD
r, sz = comm.Get_rank(), comm.Get_size()
if r == 0:
    print(f'[A] pre-barrier: size={sz}', flush=True)
comm.Barrier()
if r == 0:
    print(f'[A] barrier 1 ok', flush=True)
# exercise collective path a bit more to match Stage 3's first moments
data = comm.bcast('hello' if r == 0 else None, root=0)
comm.Barrier()
if r == 0:
    print(f'[A] bcast+barrier 2 ok', flush=True)
# one rank prints where it lives to confirm multi-node placement
print(f'[A] rank {r}/{sz} on {socket.gethostname()}', flush=True)
" 2>&1 | sed 's/^/[A] /'
rc_A=$?
t1=$(date +%s)
echo "[diag A] rc=${rc_A}  elapsed=$((t1 - t0))s"

# ----- Variant B: mpirun with MCA flags dropped (auto-detect transport)
echo
echo "=============================================================="
echo "[diag B] mpirun --oversubscribe -np N   (no MCA flags)"
echo "=============================================================="
t0=$(date +%s)
mpirun --oversubscribe -np "${N_PROCS}" \
    python -c "
from mpi4py import MPI
import socket
comm = MPI.COMM_WORLD
comm.Barrier()
if comm.Get_rank() == 0:
    print(f'[B] barrier ok  size={comm.Get_size()}', flush=True)
print(f'[B] rank {comm.Get_rank()} on {socket.gethostname()}', flush=True)
" 2>&1 | sed 's/^/[B] /'
rc_B=$?
t1=$(date +%s)
echo "[diag B] rc=${rc_B}  elapsed=$((t1 - t0))s"

# ----- Variant C: srun --mpi=pmix (SLURM-native launch)
echo
echo "=============================================================="
echo "[diag C] srun --mpi=pmix -n N   (SLURM-native)"
echo "=============================================================="
t0=$(date +%s)
srun --mpi=pmix -n "${N_PROCS}" \
    python -c "
from mpi4py import MPI
import socket
comm = MPI.COMM_WORLD
comm.Barrier()
if comm.Get_rank() == 0:
    print(f'[C] barrier ok  size={comm.Get_size()}', flush=True)
print(f'[C] rank {comm.Get_rank()} on {socket.gethostname()}', flush=True)
" 2>&1 | sed 's/^/[C] /'
rc_C=$?
t1=$(date +%s)
echo "[diag C] rc=${rc_C}  elapsed=$((t1 - t0))s"

# ----- Variant D: warmup + second mpirun (cold-start-fix candidate)
echo
echo "=============================================================="
echo "[diag D] trivial warmup mpirun, then real mpirun  (cold-start fix test)"
echo "=============================================================="
t0=$(date +%s)
mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python -c "from mpi4py import MPI; MPI.COMM_WORLD.Barrier()" \
    2>&1 | sed 's/^/[D-warm] /' || echo "[D-warm] rc=$?"

mpirun --oversubscribe -np "${N_PROCS}" \
    --mca oob_tcp_if_include ib0 \
    --mca btl_tcp_if_include ib0 \
    python -c "
from mpi4py import MPI
import socket
comm = MPI.COMM_WORLD
comm.Barrier()
if comm.Get_rank() == 0:
    print(f'[D] barrier ok  size={comm.Get_size()}', flush=True)
print(f'[D] rank {comm.Get_rank()} on {socket.gethostname()}', flush=True)
" 2>&1 | sed 's/^/[D] /'
rc_D=$?
t1=$(date +%s)
echo "[diag D] rc=${rc_D}  elapsed=$((t1 - t0))s"

echo
echo "=============================================================="
echo "[diag SUMMARY] A(prod)=${rc_A}  B(no_mca)=${rc_B}  C(srun)=${rc_C}  D(warmup+prod)=${rc_D}"
echo "=============================================================="
echo "[diag] end: $(date +%H:%M:%S.%N)"
