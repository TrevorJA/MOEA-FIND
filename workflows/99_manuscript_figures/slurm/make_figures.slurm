#!/bin/bash
#SBATCH --job-name=99_make_figures
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Stage 99 -- the SINGLE main-text figure job. Runs the per-figure
# contributor scripts (which call the src/plotting library), then
# make_figures.py assembles/re-packages everything into figures/main/.
#
#   sbatch workflows/99_manuscript_figures/slurm/make_figures.slurm [SLUG]
#
# SLUG = the stage-04 run_moea_find variant slug feeding fig05/06/07.
# Defaults to make_figures.py's own --single-site-slug default.
# Run after every contributing upstream stage has been executed.

# shellcheck source=../../_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

cd "${SLURM_SUBMIT_DIR}"

SLUG="${1:-}"
SLUG_ARG=()
[ -n "${SLUG}" ] && SLUG_ARG=(--slug "${SLUG}")

# Contributors that emit the per-variant PDFs make_figures.py repackages.
for s in dimension_sweep run_moea_find scenario_discovery_plots; do
    echo "[99/make_figures] contributor: ${s}"
    python "workflows/99_manuscript_figures/${s}.py" "${SLUG_ARG[@]}" \
        || python "workflows/99_manuscript_figures/${s}.py" \
        || echo "[99/make_figures] WARN: ${s} skipped (missing upstream?)"
done

# Final assembly (fig01..fig07 + figSI01) into figures/main/.
MF_ARG=()
[ -n "${SLUG}" ] && MF_ARG=(--single-site-slug "${SLUG}")
python workflows/99_manuscript_figures/make_figures.py "${MF_ARG[@]}"
