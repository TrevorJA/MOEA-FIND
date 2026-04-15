# Zaniolo et al. (2023) - FIND

**Full title:** FIND: A synthetic weather generator to control drought Frequency, INtensity, and Duration
**Journal:** Environmental Modelling and Software
**Authors:** Zaniolo, Fletcher, Mauter

## Summary

FIND is a MATLAB-based tool that generates synthetic streamflow traces with directly controllable and independently adjustable drought characteristics (Frequency, INtensity, Duration) while preserving observed hydrological variability and spatial correlations. Uses simulated annealing to optimize bootstrap sampling to match target FID values.

## Method

- Same simulated annealing approach as Borgomeo et al. (2015)
- Incorporates a standardized drought index for direct FID control
- Researcher specifies target FID multipliers (e.g., 1.5x frequency)
- Generator produces a single trace matching those targets

## Limitations

- **Single-objective:** Generates one scenario per run
- **No ensemble structure:** To build an ensemble, must run many times with different targets
- **No uniformity guarantee:** Manual target selection does not guarantee structured coverage
- **MATLAB only:** github.com/m-zaniolo/FIND-drought-generator

## Relevance to MOEA-FIND

MOEA-FIND generalizes FIND's single-objective FID control to a multi-objective formulation where the Pareto front provides a structured ensemble with guaranteed near-uniform coverage.
