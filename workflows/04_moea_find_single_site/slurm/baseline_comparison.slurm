#!/bin/bash
#SBATCH --job-name=moea_11_baseline_comparison
#SBATCH --output=workflows/slurm_logs/%x_%j.out
#SBATCH --error=workflows/slurm_logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
# Experiment 11: Kirsch-random-library vs MOEA-FIND baseline comparison.
# Serial; reads exp04 results.json + exp05 characteristics.npz and writes
# comparison_summary.json + fig07 scatter/convergence under
# outputs/exp11_baseline_comparison/.
#
# Defaults for input paths resolve from DEFAULT_EXPERIMENT + conventional
# locations (see workflows/04_moea_find_single_site/baseline_comparison.py).
# Override per-run via env vars:
#   sbatch --export=ALL,BORG_NFE=500000 workflows/04_moea_find_single_site/slurm/baseline_comparison.slurm

source "${SLURM_SUBMIT_DIR}/workflows/_common.sh"
moea_setup

# Resolve the exp04 variant slug to feed --moea-front. Mirrors the
# logic in 11_baseline_comparison.py::_default_moea_front but lets the
# SLURM script stay in sync when BORG_NFE is overridden.
PY_ARGS=()
if [ -n "${BORG_NFE:-}" ]; then
    PY_ARGS+=(--moea-front "$(python - <<PY
from src.experiment_config import DEFAULT_EXPERIMENT
from src.experiment_utils import make_variant_slug
from pathlib import Path
import os
cfg = DEFAULT_EXPERIMENT
nfe = int(os.environ.get("BORG_NFE", cfg.nfe))
slug = make_variant_slug(mode=cfg.dv_mode, n_years=cfg.n_years_out, nfe=nfe,
                         seed=cfg.seed, constrained=cfg.constraints_json is not None)
print(Path("outputs/exp04_kirsch_single_site") / slug / "results.json")
PY
)")
fi

python workflows/04_moea_find_single_site/baseline_comparison.py "${PY_ARGS[@]}"
echo "exit: $?"
