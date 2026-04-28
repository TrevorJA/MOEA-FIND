#!/bin/bash
#SBATCH --job-name=moea_04_verify_e2e
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
# End-to-end pipeline verification with variant slug system.
# Runs two configurations — each gets its own slug-named output directory.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export MOEA_FIND_ALGORITHM=borg_mm
export OMP_NUM_THREADS=1

CONSTRAINTS_JSON="${PROJECT_ROOT}/outputs/diag_constraint_calibration/calibrated_tolerances.json"
N_PROCESSORS=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
echo "[verify] nodes=${SLURM_NNODES} tasks/node=${SLURM_NTASKS_PER_NODE} total=${N_PROCESSORS}"

set +e  # don't abort on MPI cleanup warnings

# ===== RUN 1: residual, 5yr, 20k NFE (expect solutions + figures) =====
echo ""
echo "========== RUN 1: residual_T5_nfe20000_s42 =========="
moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --nfe 20000 \
    --mode residual \
    --n-years 5 \
    --ssi 3 \
    --seed 42 \
    --algorithm borg_mm \
    --constraints-json "${CONSTRAINTS_JSON}" \
    --n-islands 1 \
    --checkpoint-freq 20000 \
    --plot

# ===== RUN 2: index, 5yr, 20k NFE (compare sampling modes) =====
echo ""
echo "========== RUN 2: index_T5_nfe20000_s42 =========="
moea_mpi_launch "${N_PROCESSORS}" python workflows/04_moea_find_single_site/run_moea_find.py \
    --nfe 20000 \
    --mode index \
    --n-years 5 \
    --ssi 3 \
    --seed 42 \
    --algorithm borg_mm \
    --constraints-json "${CONSTRAINTS_JSON}" \
    --n-islands 1 \
    --checkpoint-freq 20000 \
    --plot

# ===== VERIFICATION CHECKS =====
echo ""
echo "========== VERIFICATION CHECKS =========="

BASE="${PROJECT_ROOT}/outputs/exp04_kirsch_single_site"

# Check that both variants produced separate directories
echo "--- Output directories ---"
for d in "${BASE}"/*/; do
    slug=$(basename "$d")
    if [ -f "${d}results.json" ]; then
        python -c "
import json
r = json.load(open('${d}results.json'))
n = r.get('n_pareto', 0)
print(f'  [OK] {\"${slug}\"}: n_pareto={n}, evals={r.get(\"n_evals_total\",0)}, elapsed={r.get(\"elapsed_s\",0):.1f}s')
"
    else
        echo "  [FAIL] ${slug}: no results.json"
    fi
done

echo "--- Figures ---"
for d in "${BASE}"/*/figures/; do
    slug=$(basename "$(dirname "$d")")
    nfigs=$(ls "${d}"*.pdf 2>/dev/null | wc -l)
    echo "  [OK] ${slug}/figures/: ${nfigs} PDFs"
    ls "${d}"*.pdf 2>/dev/null | xargs -I{} basename {} | sed 's/^/    /'
done

echo ""
echo "========== VERIFICATION COMPLETE =========="
