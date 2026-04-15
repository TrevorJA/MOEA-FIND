# Parametric Streamflow Generation: Beyond Bootstrap Resampling

## The Problem

Pure block bootstrap (Kirsch-Nowak) can only recombine observed historical patterns. Generated droughts cannot exceed what was observed. For stress testing, we need methods that extend beyond the historical envelope.

## State-of-the-Art Approaches

### Kappa Distribution for Monthly Flow Marginals

Svensson et al. (2017, WRR) established the four-parameter kappa distribution as a flexible choice for monthly streamflow. Bounded below at zero, with flexible tail behavior for extremes. Widely used for drought indicator applications.

**Reference:** Svensson, C., et al. (2017). Statistical distributions for monthly aggregations of precipitation and streamflow in drought indicator applications. WRR. doi:10.1002/2016WR019276

### Vine Copulas for Temporal and Spatial Dependence

When sampling from parametric marginals (rather than resampling historical values), temporal and spatial correlation must be imposed separately. Vine copulas are the current standard:

- **D-vine copulas** with lag-1 and lag-2 terms capture temporal autocorrelation
- **Periodic vine copulas** account for seasonally varying dependence
- **PAR(p)-vine models** combine periodic autoregression with vine copulas

**Key references:**
- Brechmann, E.C., et al. (2017). PAR(p)-vine copula based model for stochastic streamflow scenario generation. Stoch. Environ. Res. Risk Assess.
- Li, C., et al. (2023). Mixed D-vine copula-based conditional quantile model for stochastic monthly streamflow simulation. J. Hydro-environ. Res.
- Periodic vine copula autoregressive model (2019). Water Resour. Manage.

### Kernel Density Estimation (KDE)

Nonparametric smoothing of empirical distributions. Extends support beyond observed range via kernel bandwidth while staying close to data.

- Mixed copula-KDE framework (2025, Stoch. Environ. Res. Risk Assess.) for joint streamflow modeling

### Phase Randomization (Papalexiou and Koutsoyiannis)

Decomposes streamflow into amplitude and phase in Fourier space. Randomizes phase, re-scales using kappa distribution parameters. Preserves spectral properties while allowing parametric tail extension.

**Reference:** Papalexiou, S.M. (2019). Stochastic simulation of streamflow time series using phase randomization. HESS.

### Bootstrapped CDF Uncertainty

Bootstrap the historical data to generate ensemble of CDF estimates with confidence bands. Sample from the ensemble rather than a single point estimate.

**Reference:** BEUM approach (2022, Stoch. Environ. Res. Risk Assess.)

## Relevance to MOEA-FIND

For the publication, a two-track approach:
1. **Bootstrap (Kirsch-Nowak):** Baseline formulation. DVs are sample indices. Limited to historical envelope.
2. **Parametric (kappa + vine copula):** Advanced formulation. DVs are CDF probabilities. Extends beyond historical.

The comparison between tracks is itself a contribution: does parametric sampling meaningfully expand the achievable drought characteristic space?
