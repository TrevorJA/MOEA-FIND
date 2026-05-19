#!/bin/bash
#SBATCH --job-name=99_make_si_figures
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Stage 99-SI -- the SINGLE supporting-information figure job. Iterates
# every script in workflows/99_supporting_info_figures/ (each calls the
# src/plotting library). Best-effort: a script whose upstream outputs
# are absent is logged and skipped, not fatal, so one missing stage
# never blocks the rest.
#
#   sbatch workflows/99_supporting_info_figures/slurm/make_si_figures.slurm [SLUG]
#
# SLUG = upstream variant slug for scripts that take --slug. Scripts with
# a different arg signature fall back to a no-arg call. Per-script slugs
# that differ from the default can be refined here as the SI matures.

# shellcheck source=../../_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

cd "${SLURM_SUBMIT_DIR}"

SLUG="${1:-}"
DIR=workflows/99_supporting_info_figures
fail=0; ok=0; skip=0

for f in "${DIR}"/*.py; do
    [ "$(basename "${f}")" = "__init__.py" ] && continue
    echo "[99-SI] === ${f} ==="
    if [ -n "${SLUG}" ] && python "${f}" --slug "${SLUG}"; then
        ok=$((ok+1))
    elif python "${f}"; then
        ok=$((ok+1))
    else
        echo "[99-SI] SKIP/FAIL ${f} (missing upstream or needs explicit args)"
        skip=$((skip+1))
    fi
done

echo "[99-SI] done: ${ok} produced, ${skip} skipped/failed (non-fatal)"
