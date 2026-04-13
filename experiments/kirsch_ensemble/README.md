# Hydrology-Based Synthetic Ensemble Experiments

This directory contains experiments that generate synthetic streamflow traces using data-driven methods (Kirsch bootstrap and parametric generators) coupled with multi-objective optimization.

## Scripts

### `run_kirsch_experiment.py`

**Kirsch nonparametric bootstrap with Borg wrapper**

Tests the `KirschBorgWrapper` class with both "index" and "residual" decision variable modes using SSI-based drought objectives. The Kirsch method preserves temporal autocorrelation, cross-year correlation (Dec-Jan), and seasonal structure via Cholesky decomposition and normal-score transforms.

This experiment validates that the wrapper correctly maps Borg DVs to indices/residuals and generates physically plausible synthetic traces.

Usage:
```bash
python experiments/kirsch_ensemble/run_kirsch_experiment.py --nfe 5000
python experiments/kirsch_ensemble/run_kirsch_experiment.py --mode residual --nfe 5000
```

### `run_parametric_experiment.py`

**Parametric CDF + vine copula generator**

Tests the parametric generator (kappa4 marginals + D-vine copula) for drought characteristic space exploration. Uses threshold-based drought metrics (frequency, duration, magnitude, severity).

Operates on USGS gauge 01423000 (West Branch Delaware at Walton, Cannonsville inflow).

Usage:
```bash
python experiments/kirsch_ensemble/run_parametric_experiment.py --nfe 5000
```

## Output

Results are saved to `../../outputs/kirsch_ensemble/` and `../../outputs/parametric/` respectively, including:
- JSON files with Pareto front solutions
- Figure directories with diagnostic plots
