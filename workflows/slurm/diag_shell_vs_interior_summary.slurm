#!/bin/bash
#SBATCH --job-name=moea_diag_shell_vs_interior_summary
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
# Post-process the per-K diag_shell_vs_interior runs into the
# dimension-sweep summary figure (figSI_shell_interior_sweep.pdf) and
# markdown table (outputs/diag_shell_vs_interior/sweep_table.md).
#
# Intended to be submitted with --dependency=afterok:<array_job_id>
# immediately after the array job in diag_shell_vs_interior.slurm.
# Safe to run standalone any time the per-K results.json files exist.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

echo "[summary.slurm] start: $(date)"
python workflows/diagnostics/diag_shell_vs_interior_summary.py
echo "[summary.slurm] end:   $(date)"
