#!/bin/bash
#SBATCH --job-name=04_first_event_ssi3_t10_fixed
#SBATCH --output=workflows/slurm_logs/%x_%A.out
#SBATCH --error=workflows/slurm_logs/%x_%A.err
#SBATCH --time=08:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=30
#SBATCH --exclusive
# Stage 04 / T=10 first-SSI3-drought, fixed-anti-ideal 5-objective set.
# Single 200k-NFE MM Borg run on 8 exclusive nodes x 30 ranks = 240 ranks.
# Objectives (preset first_event_ssi3_t10_fixed), first critical SSI-3
# drought only, NO duration:
#   severity (D*=6.75, eps=0.1) | magnitude (D*=120, eps=1)
#   | peak_month (D*=13, eps=1) | onset_rate (D*=6.75, eps=0.1)
#   | rising_limb (D*=1.1, eps=0.1)
# severity/magnitude/onset = user tight targets (4.5/80/4.5) x 1.5 headroom.
# All AntiIdealRule.CONSTANT — D* fixed, not historical/cyclic headroom.
# Hard ">=1 critical SSI-3 event" constraint auto-gated by first_event_
# name prefix; DV-uniformity Anderson-Darling constraint as production.
# Calibration prerequisite:
#   outputs/02_calibration/dv_uniformity_calibration/residual/calibrated_dv_tolerances.json
# Runs alongside first_event_ssi3_t10_iter1.slurm (no clobber: --slug-suffix
# fixed vs iter1 gives distinct output dirs at the same mode/T/nfe/seed).

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[04/first_event_ssi3_t10_fixed] cores=${N_PROCESSORS}"

moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --config workflows/04_moea_find_single_site/configs/first_event_ssi3_t10_fixed.yaml \
    --slug-suffix fixed
