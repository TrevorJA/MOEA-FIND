# Analytic Proof of Concept Experiments

This directory contains pure analytic tests of the MOEA-FIND methodology on synthetic decision spaces. These tests do not involve hydrology or streamflow data.

## Active Scripts

### `run_analytic_poc.py`

Runs 2D and 3D analytic test cases:

- **Experiment 1.1: 2D Analytic Test (Replication)** - Replicate the original 2D proof-of-concept from the proposal:
  - J_1 = X_1, J_2 = X_2, J_3 = ||X - X*||_1
  - X in [-3, 3]^2, X* = (3, 3)
  - Verify near-uniform coverage of (X_1, X_2) from Borg's Pareto front
  - Compute discrepancy metrics and compare to LHS baseline

- **Experiment 1.2: 3D Analytic Extension** - Extend to 3 objectives + Manhattan norm (4 total):
  - J_1 = X_1, J_2 = X_2, J_3 = X_3, J_4 = ||X - X*||_1
  - Verify uniformity on the 2-simplex
  - Key test of generalizability

- **Experiment 1.3: Epsilon Sensitivity** - Vary epsilon values and measure coverage quality.

### `plot_poc_results.py`

Generate diagnostic figures from analytic POC results.

## Hydrology Experiments (Moved)

Hydrology-based synthetic ensemble experiments have been moved to `../kirsch_ensemble/`:
- `run_kirsch_experiment.py` (previously `run_kirsch_poc.py`)
- `run_parametric_experiment.py` (previously `run_parametric_poc.py`)

Results are in `../../outputs/kirsch_ensemble/` and `../../outputs/parametric/`.

## Deprecated Scripts

The following scripts are deprecated and have been superseded:
- `run_bootstrap_poc.py` — old bootstrap experiment (removed)
- `plot_bootstrap_diagnostics.py` — old bootstrap plotting (removed)
- `run_ssi_poc.py` — old SSI test (removed)
- `run_kirsch_poc.py` — moved to `../kirsch_ensemble/run_kirsch_experiment.py`
- `run_parametric_poc.py` — moved to `../kirsch_ensemble/run_parametric_experiment.py`
