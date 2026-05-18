#!/bin/bash
#SBATCH --job-name=04_first_event_ssi3_t10_smoke
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks=120
#SBATCH --mem-per-cpu=2G
#SBATCH --exclusive
# Stage 04 / T=10 first-SSI3-drought — Phase D smoke wiring test.
# Verifies:
#   * preset 'first_event_ssi3_t10' resolves to 5 first-event metrics
#   * at-least-one-event hard constraint fires for empty-event candidates
#   * DV-uniformity AD constraint loads cannonsville_T10_ad tolerance
#   * Borg runs to completion and writes results.json
# Production 120-rank allocation per the project parallelize-by-default
# rule; _auto_islands picks 7 islands at this scale. NFE kept low for
# smoke; the iter1 production script will rerun with 200k NFE.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=${SLURM_NTASKS}
echo "[04/first_event_ssi3_t10_smoke] tasks=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/first_event_ssi3_t10.yaml \
    --nfe 20000 \
    --slug-suffix smoke
