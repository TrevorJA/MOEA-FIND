# DRB Case Study

Application of MOEA-FIND to the Delaware River Basin for the NYCOptimization policy re-evaluation.

## Setup

- Multi-site generation: Cannonsville, Pepacton, Neversink inflows + lateral flows
- Drought characteristics: frequency, duration, intensity (3 objectives + Manhattan = 4 total)
- Block bootstrap with shared indices for spatial correlation
- Generated ensemble used to re-evaluate Pareto-optimal policies from NYCOptimization

## Dependencies

- SynHydro (for historical data loading and preprocessing)
- Pywr-DRB (for policy re-evaluation simulation)
- Borg MOEA (for optimization)
