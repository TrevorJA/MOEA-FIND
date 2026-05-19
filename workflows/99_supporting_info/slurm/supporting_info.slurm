#!/bin/bash
#SBATCH --job-name=99_supporting_info
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Stage 99-SI -- the SINGLE supporting-information COMPUTE job. Runs the
# SI-only diagnostics relocated out of 02_calibration during the
# spring-cleaning reorg (kirsch convergence + Kirsch-wrapper geometry/
# fidelity). Each is compute-only; figures come from the matching
# scripts in workflows/99_supporting_info_figures/.
#
# The --library / --moea-front paths below are run-specific; update the
# slug to match the produced stage-03 / stage-04 outputs before submit.

# shellcheck source=../../_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1

LIBRARY="${1:-outputs/03_kirsch_library/build_library/library__N=10000__T=20__s=42/characteristics.npz}"
MOEA_FRONT="${2:-outputs/04_moea_find_single_site/run_moea_find/$(ls outputs/04_moea_find_single_site/run_moea_find 2>/dev/null | head -1)/results.json}"

echo "[99-SI] kirsch_convergence"
python workflows/99_supporting_info/kirsch_convergence.py \
    --library "${LIBRARY}" --moea-front "${MOEA_FRONT}" \
    --n-seeds 10 --seed 42 \
    || echo "[99-SI] WARN: kirsch_convergence skipped (check --library/--moea-front)"

echo "[99-SI] wrapper_geometry"
python workflows/99_supporting_info/wrapper_geometry.py \
    --seed 42 --n-sweep-points 400 --n-years 20

echo "[99-SI] wrapper_fidelity"
python workflows/99_supporting_info/wrapper_fidelity.py \
    --n-traces 1000 --n-years 20 --ssi 3 --seed 42
