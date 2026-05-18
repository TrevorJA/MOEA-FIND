#!/bin/bash
#SBATCH --job-name=02_figures_abc
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Stage 02 / figures A, B, C — manuscript figures supporting the joint
# (K*, T*) decision. Reads Stage 1+2+3 outputs and (optionally) the
# Stage 4 MOEA-FIND archive via PARETO_ARCHIVE env var.
#
# Submit:
#   sbatch workflows/02_calibration/slurm/plots/figures_a_b_c.slurm
#   sbatch --export=ALL,PARETO_ARCHIVE=outputs/04_moea_find_single_site/<slug>/results.json \
#          workflows/02_calibration/slurm/plots/figures_a_b_c.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

PARETO_ARCHIVE="${PARETO_ARCHIVE:-}"
ARCHIVE_FLAG=""
if [[ -n "${PARETO_ARCHIVE}" ]]; then
    ARCHIVE_FLAG="--pareto-archive ${PARETO_ARCHIVE}"
    echo "[figures_abc] using Pareto archive: ${PARETO_ARCHIVE}"
fi

echo "[figures_abc] Figure A — metric stability across T"
python src/plotting/02_calibration/figure_a_metric_stability.py \
    --T-grid 5 10 20 30 ${ARCHIVE_FLAG}

echo "[figures_abc] Figure B — correlation clustermaps"
python src/plotting/02_calibration/figure_b_corr_clustermaps.py \
    ${ARCHIVE_FLAG}

echo "[figures_abc] Figure C — drought-space coverage"
python src/plotting/02_calibration/figure_c_drought_space.py \
    ${ARCHIVE_FLAG}
