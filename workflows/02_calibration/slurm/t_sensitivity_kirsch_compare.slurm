#!/bin/bash
#SBATCH --job-name=02_t_sensitivity_kirsch_compare
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# Stage 02 / t_sensitivity_kirsch_compare -- Kirsch vs historical fidelity.
# Run AFTER both Stage 1 historical array AND Stage 2 build_library_extended
# array complete. Submit with:
#   sbatch --dependency=afterok:<hist_jid>:<kirsch_jid> \
#          workflows/02_calibration/slurm/t_sensitivity_kirsch_compare.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

echo "[02/t_sensitivity_kirsch_compare] running"
python workflows/02_calibration/t_sensitivity_kirsch_compare.py \
    --T-grid 5 10 20 30 \
    --n-traces 10000 \
    --seed 42
