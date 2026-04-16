#!/bin/bash
#SBATCH --job-name=moea_09_drb_policy_reeval
#SBATCH --output=workflows/slurm/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm/slurm_logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
# Manuscript: §7.3 Policy Re-evaluation (Figure 9)
# Experiment:  09 — MOEA-FIND ensemble re-evaluated on NYCOptimization policies
# Parallelism: serial Pywr-DRB driver; parallelism inside Pywr-DRB if supported
# STATUS:      SCAFFOLDED; pending script 08 ensemble and NYCOptimization handoff.

# shellcheck source=./_common.sh
source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

python workflows/experiments/09_drb_policy_reeval.py \
    --ensemble outputs/exp08_drb_multisite/ensemble.npz \
    --policies "${NYC_POLICIES_PATH:-../NYCOptimization/outputs/pareto_policies.json}"
