# Wheeler et al. (2025) — Multisite Nonparametric Stochastic Streamflow Generation

**Full title:** Multisite Nonparametric Stochastic Streamflow Generation for the Eastern Nile Basin
**Journal:** Journal of Hydrologic Engineering (ASCE)
**Authors:** Wheeler, Simpson, Borgomeo, Hall

## Summary

Extends the Borgomeo/Zaniolo simulated annealing approach to multiple sites. Key differences from Zaniolo (2024):
- Focus on **multiple sites** using cross-site correlation matrix in the objective function
- **Hurst coefficient** included in the objective function (long-range dependence)
- Same simulated annealing optimizer

## Relevance to MOEA-FIND

Shows that search-based generation can handle multi-site requirements. The cross-site correlation can be handled either as an objective or a constraint. MOEA-FIND's block bootstrap with shared indices naturally preserves spatial correlation (a simpler approach than including it in the objective).
