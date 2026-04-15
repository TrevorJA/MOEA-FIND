# Borgomeo et al. (2015) - Risk-based water resources planning

**Full title:** Risk-based water resources planning: Incorporating probabilistic nonstationary climate uncertainties
**Journal:** Water Resources Research
**Key contribution:** First to frame streamflow generation as a combinatorial optimization problem

## Summary

Borgomeo et al. proposed treating streamflow generation as an optimization problem where the objective is to match one or more user-specified statistics of the resulting flows. Uses a stochastic search algorithm (simulated annealing) to find bootstrap resamplings of historical data that produce traces with desired properties.

## Relevance to MOEA-FIND

This is the most direct intellectual predecessor. The key limitation is that it is **single-objective**: generates one trace per optimization run to match one target. MOEA-FIND extends this to multi-objective, where the Pareto front is the ensemble.

## Key Quote (from proposal slide 4)

"Borgomeo et al. (2015) provided perhaps the most general approach, in which streamflow generation is treated as a combinatorial optimization problem to match one or more user specified statistics of the resulting flows."
