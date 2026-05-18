#!/bin/bash
# Orchestrator blueprint for the joint metric-and-T sensitivity workflow.
#
# Submits Stages 1, 2, 3 with SLURM dependencies so each stage waits for
# its inputs. NOT auto-run on the login node — read it, copy/paste the
# `sbatch` lines, or `bash workflows/02_calibration/run_t_sensitivity.sh`
# from the project root.
#
# Stage 4 (MOEA-FIND validation) is intentionally NOT submitted here:
# the user reviews `outputs/02_calibration/decision_matrix/pareto_front_KxT.json`
# before committing 120-rank Borg compute. Same for Stage 5 figures.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

# Stage 1 — historical-block per-T screening.
HIST_JID=$(sbatch --parsable \
    workflows/02_calibration/slurm/t_sensitivity_historical.slurm)
echo "Stage 1 historical: ${HIST_JID}"

AGG_JID=$(sbatch --parsable --dependency=afterok:"${HIST_JID}" \
    workflows/02_calibration/slurm/t_sensitivity_aggregate.slurm)
echo "Stage 1 aggregate:  ${AGG_JID}"

# Stage 2 — Kirsch ensemble + fidelity comparison.
KIRSCH_JID=$(sbatch --parsable \
    workflows/03_kirsch_library/slurm/build_library_extended.slurm)
echo "Stage 2 Kirsch:     ${KIRSCH_JID}"

KCOMP_JID=$(sbatch --parsable \
    --dependency=afterok:"${KIRSCH_JID}":"${AGG_JID}" \
    workflows/02_calibration/slurm/t_sensitivity_kirsch_compare.slurm)
echo "Stage 2 compare:    ${KCOMP_JID}"

# Stage 3 — eval cost timing + decision matrix.
TIMING_JID=$(sbatch --parsable \
    workflows/02_calibration/slurm/eval_cost_timing.slurm)
echo "Stage 3 timing:     ${TIMING_JID}"

DECISION_JID=$(sbatch --parsable \
    --dependency=afterok:"${KCOMP_JID}":"${TIMING_JID}" \
    workflows/02_calibration/slurm/decision_matrix.slurm)
echo "Stage 3 decision:   ${DECISION_JID}"

cat <<EOF

---
All stages submitted. Once ${DECISION_JID} completes, review:
  outputs/02_calibration/decision_matrix/pareto_front_KxT.json
  outputs/02_calibration/decision_matrix/decision_matrix.csv

Then register the recommended K* preset in src/metrics/drought_metrics.py
and submit Stage 4 + Stage 5 (figures) manually.
EOF
