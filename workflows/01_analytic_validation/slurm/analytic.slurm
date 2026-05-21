#!/bin/bash
#SBATCH --job-name=01_analytic
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
# Stage 01 / analytic validation -- MM Borg via MPI on the K-dimensional
# analytic benchmark. Array index = K - 2, so tasks 0..4 run K=2..6.
# Earlier per-K analytic_2d/3d drivers were retired -- the K=2 and K=3
# slices of this sweep produce identical archives at higher NFE, and
# the SI-1 hyperplane + Fig 4 dimension figures both read from
# outputs/01_analytic_validation/dimension_sweep/k{K}/.
#
# The 36-cell epsilon x NFE sweep stays separate (eps_nfe_sweep.slurm,
# different rank count); aggregate via eps_nfe_aggregate.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

K=$((SLURM_ARRAY_TASK_ID + 2))   # tasks 0..4 -> K=2..6

# Per-K (eps, NFE) -- low-K runs get more NFE so the K=2/K=3 slices
# match the resolution the dedicated analytic_2d/3d drivers produced.
case "${K}" in
  2) EPS=0.10; NFE=100000 ;;
  3) EPS=0.15; NFE=50000  ;;
  4) EPS=0.20; NFE=10000  ;;
  5) EPS=0.25; NFE=10000  ;;
  6) EPS=0.30; NFE=10000  ;;
esac

CMD=(python workflows/01_analytic_validation/dimension_sweep.py
     --k "${K}" --nfe "${NFE}" --epsilon "${EPS}"
     --feasible-shape cube --seed 42)

echo "[01/analytic] K=${K} NFE=${NFE} EPS=${EPS} (ntasks=${SLURM_NTASKS})"
moea_mpi_launch "${SLURM_NTASKS}" "${CMD[@]}"
