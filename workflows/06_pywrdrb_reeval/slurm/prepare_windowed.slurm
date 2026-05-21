#!/bin/bash
#SBATCH --job-name=06_prepare_windowed
#SBATCH --output=workflows/slurm_logs/%x_%A_%a.out
#SBATCH --error=workflows/slurm_logs/%x_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --array=0-9
# Stage 06 (WINDOWED) per-scenario INPUT PREPARATION for the COARSE
# stage-04 archive (DV -> daily multisite HDF5 + per-scenario window).
#
# Split out of policy_reeval_windowed.slurm (2026-05-20) so the heavy
# Kirsch-replay + Nowak-disagg + KDE chain runs once per archive and
# the simulation step can be re-run without redoing it. Idempotent:
# scenarios with status_prep.json ok=True are skipped on re-submit
# (pass --overwrite to force re-prep).
#
# Smart hybrid parallelism (matches policy_reeval_windowed.slurm so the
# two stages share allocation profile):
#   * SLURM array (10 tasks) chunks the ~1009 scenarios
#   * Within each array task: mpirun -np 30, one core per rank
#   * intra-node shared-memory only (no inter-node ORTE)
# 10 tasks x 30 ranks = 300 workers over ~1009 scenarios ~= 3-4 each.

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Production first-event slug emitted by run_moea_find.py with
# configs/first_event_ssi3_t10.yaml. Override by passing the full
# results.json path as $1, or by editing this default in-place.
MOEA_SLUG="${MOEA_SLUG:-moea__mode=residual__T=10__nfe=200k__ssi=3__metrics=first_event_ssi3_t10__cons=dv-l2__st=ad__s=42}"
PARETO_RESULTS="${1:-outputs/04_moea_find_single_site/run_moea_find/${MOEA_SLUG}/results.json}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

N_PROCS="${SLURM_NTASKS:-30}"
echo "[06prep.slurm] array task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}, ${N_PROCS} ranks (1 core/rank)"
echo "[06prep.slurm] Pareto: ${PARETO_RESULTS}"
echo "[06prep.slurm] Start: $(date)"

mpirun --oversubscribe -np "${N_PROCS}" \
    --mca pml ob1 \
    --mca btl self,vader \
    python "${SLURM_SUBMIT_DIR}/workflows/06_pywrdrb_reeval/prepare_windowed.py" \
        --pareto-results "${PARETO_RESULTS}" \
        --site-label cannonsville \
        --ssi 3 \
        --demand-source constant_max \
        --use-mpi

echo "[06prep.slurm] End: $(date)"
