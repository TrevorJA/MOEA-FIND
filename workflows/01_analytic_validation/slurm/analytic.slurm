#!/bin/bash
#SBATCH --job-name=01_analytic
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --array=0-6
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Stage 01 / analytic validation -- combined launcher (MM Borg via MPI).
# Replaces the former analytic_2d / analytic_3d / dimension_sweep slurms.
# Array dispatch:
#   task 0       analytic_2d   (Fig SI-1)
#   task 1       analytic_3d   (Fig SI-1)
#   task 2..6    dimension_sweep K=2..6  (Fig 4)
# The 36-cell epsilon x NFE sweep stays separate (eps_nfe_sweep.slurm,
# different rank count); aggregate via eps_nfe_aggregate.slurm.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

T="${SLURM_ARRAY_TASK_ID:-0}"
case "${T}" in
  0) CMD=(python workflows/01_analytic_validation/analytic_2d.py
          --nfe 100000 --seed 42 --epsilon 0.06) ;;
  1) CMD=(python workflows/01_analytic_validation/analytic_3d.py
          --nfe 50000 --seed 42 --epsilon 0.15) ;;
  *) K=$((T))                       # tasks 2..6 -> K=2..6
     case "${K}" in
       2) EPS=0.10 ;; 3) EPS=0.15 ;; 4) EPS=0.20 ;;
       5) EPS=0.25 ;; *) EPS=0.30 ;;
     esac
     CMD=(python workflows/01_analytic_validation/dimension_sweep.py
          --k "${K}" --nfe 5000 --epsilon "${EPS}"
          --feasible-shape cube --seed 42) ;;
esac

echo "[01/analytic] task=${T} -> ${CMD[*]} (ntasks=${SLURM_NTASKS})"
moea_mpi_launch "${SLURM_NTASKS}" "${CMD[@]}"
