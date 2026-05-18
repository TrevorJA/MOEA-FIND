#!/bin/bash
#SBATCH --job-name=02_recompute_bounded_candidates
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / recompute_bounded_candidates -- DD-15c Phase 2.
# Recomputes the 48 bounded candidate metrics on the historical record
# and on the four "before" Pareto archives (residual+AD, +AD-iter2, +KS,
# index+AD). No MPI, no Borg. Outputs feed the Phase 3 K-set selection
# and the upgraded calibration figures.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/recompute_bounded_candidates] running"
python workflows/02_calibration/recompute_bounded_candidates.py
