#!/bin/bash
#SBATCH --job-name=moea_diag_shell_vs_interior
#SBATCH --output=workflows/slurm/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%A_%a.err
#SBATCH --array=2-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
# Analytic K-dim interior-coverage sanity check for the DD-11 L1 Device.
#
# Submits as a SLURM **job array** so the K=2..6 runs execute in
# parallel (each is ~15 s on 1 CPU). A separate afterok summary job
# rolls the per-K results.json into the sweep figure once all array
# tasks complete.
#
# Runs workflows/diagnostics/diag_shell_vs_interior.py on a K-dim cube
# feasibility inside the bounding box [-3, 3]^K (override via
# DIAG_SHAPE=ball). Uses src.objectives.analytic_objectives — the same
# DD-11 formulation as drought_objectives.
#
# Per-K outputs at outputs/diag_shell_vs_interior/k{K}/; figures at
# figures/figSI_shell_interior_k{K}.pdf; sweep summary at
# figures/figSI_shell_interior_sweep.pdf.
#
# Submit (array + summary together):
#   ARR_ID=$(sbatch --parsable workflows/slurm/diag_shell_vs_interior.slurm)
#   sbatch --dependency=afterok:${ARR_ID} workflows/slurm/diag_shell_vs_interior_summary.slurm
#
# Override NFE / epsilon / shape / seed / K-list from sbatch via env vars:
#   sbatch --export=ALL,DIAG_NFE=10000,DIAG_EPS=0.20 workflows/slurm/diag_shell_vs_interior.slurm
#   sbatch --export=ALL,DIAG_SHAPE=ball workflows/slurm/diag_shell_vs_interior.slurm
#   sbatch --array=3-5 workflows/slurm/diag_shell_vs_interior.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

export OMP_NUM_THREADS=1

: "${DIAG_NFE:=5000}"
: "${DIAG_EPS:=0.25}"
: "${DIAG_SHAPE:=cube}"
: "${DIAG_SEED:=42}"

K="${SLURM_ARRAY_TASK_ID}"

echo "[diag.slurm] K=${K} NFE=${DIAG_NFE} eps=${DIAG_EPS} shape=${DIAG_SHAPE} seed=${DIAG_SEED}"
echo "[diag.slurm] start: $(date)"

python workflows/diagnostics/diag_shell_vs_interior.py \
    --k "${K}" \
    --nfe "${DIAG_NFE}" \
    --epsilon "${DIAG_EPS}" \
    --feasible-shape "${DIAG_SHAPE}" \
    --seed "${DIAG_SEED}"

echo "[diag.slurm] end: $(date)"
