# MOEA-FIND

**Multi-Objective Evolutionary Algorithm for Structured Synthetic Drought Generation**

A novel method that couples Borg MOEA with block bootstrap streamflow generation (Kirsch-Nowak) to produce ensembles of synthetic traces with structured, near-uniform coverage of drought characteristic space. Extends the single-objective search-based generation of Borgomeo et al. (2015) and Zaniolo et al. (2023, FIND) to a multi-objective formulation where the Pareto front IS the drought ensemble.

---

## Current State (2026-04-13, Session 14 Complete)

**Key accomplishments:**
- **Fabricated literature corrected:** Removed unverified claims about Quinn et al., Hadjimichael et al., and Herman et al. using post-hoc drought space subsampling. Verified via Zotero and web search: all cited literature operates in parameter or input space, not drought outcome space. No prior publication uses space-filling subsampling in drought characteristic space.
- **Bonham et al. (2024) added** as the closest published precedent (subsamples in uncertainty input space for Colorado River Basin, R2=0.77-0.91 for ranking preservation; not in drought outcome space).
- **3D analytic validation COMPLETE:** 1362 Pareto solutions, hyperplane verified to machine precision (std~10^-16). Coverage metrics: NN_CV=0.42 (Pareto) vs 0.37 (LHS) vs 0.28 (Sobol); L2*=0.038 (Pareto) vs 0.006 (LHS) vs 0.001 (Sobol). See DD-10 for analysis: epsilon-dominance produces less uniform coverage than QMC in unconstrained analytic test, but real advantage emerges in constrained feasible drought regions (hydrology).
- **6 publication-quality figures GENERATED:** fig1 (Manhattan norm concept), fig2 (3D simplex coverage), fig3 (coverage comparison: Pareto vs LHS vs Sobol vs Random), fig4 (hyperplane verification), fig5 (epsilon sensitivity), fig6 (NN distance distributions). PNG + PDF formats ready.
- **MM Borg MOEA resources identified:** Feb 2025 Water Programming posts on serial Python wrapper; Aug 2025 post on MM Borg + checkpointing + MOEAFramework 5.0; passNFE_ALH_PyCheckpoint branch confirmed.
- **Library-and-subsample baseline implemented:** `src/library.py` and `experiments/kirsch_ensemble/run_library_baseline.py` ready for 10K+ trace generation and LHS/Sobol subsampling in drought characteristic space.
- **Design decisions finalized:** DD-06 (KirschBorgWrapper via SynHydro), DD-03 (Kirsch B=1 + Cholesky pipeline), DD-09/DD-10 (coverage baselines and epsilon-dominance analysis).

**Next immediate steps:**
1. Run Phase 2 Experiment 2.1 (single-site Kirsch-based MOEA-FIND with serial Borg or EpsNSGAII, SSI-3 drought objectives: duration + avg_severity + Manhattan norm)
2. Run Experiment 2.3 (coverage comparison: MOEA-FIND vs library+LHS subsample vs FIND-with-grid vs Sobol baseline vs random)
3. Generate 10K+ Kirsch library on HPC or local high-memory machine for baseline comparison
4. Transition to MM Borg for parallel hydrology experiments (HPC access required; Python wrapper setup per Feb 2025 Water Programming posts)

## Project Structure

```
MOEA-FIND/
├── src/                    # Core library
│   ├── __init__.py
│   ├── objectives.py       # Drought metrics, SSI, Manhattan norm
│   ├── kirsch_wrapper.py   # KirschBorgWrapper (wraps SynHydro)
│   ├── library.py          # Library-and-subsample baseline
│   ├── parametric.py       # Kappa4 + D-vine copula generator
│   ├── constraints.py      # Plausibility constraints
│   ├── analysis.py         # Coverage metrics, space-filling baselines
│   ├── data.py             # USGS data loading utilities
│   └── plotting/           # Modular diagnostic plotting
├── experiments/            # Reproducible experiment scripts
│   ├── proof_of_concept/   # Analytic and single-site POC
│   ├── kirsch_ensemble/    # Kirsch-based MOEA-FIND experiments
│   └── drb_case_study/     # DRB multi-site application (planned)
├── notes/                  # Research notes, design decisions
│   └── literature/         # Paper summaries
└── outputs/                # Generated data (gitignored)
```

## Key Concepts

- **Decision variables:** Monthly bootstrap sample indices (or CDF probabilities), array of size N_years * 12
- **Objectives:** k drought characteristics (frequency, duration, intensity, ...) + L1 Manhattan norm from anti-ideal
- **Manhattan norm trick:** J_{k+1} = ||D_i - D*_i||_1 forces the Pareto front onto a hyperplane. Borg's epsilon-dominance tiles this hyperplane near-uniformly. Projection to k-dim drought space yields structured coverage.
- **Constraints:** Plausibility of generated traces (autocorrelation, non-drought period statistics)
- **Output:** Pareto front = ensemble of synthetic streamflow traces with near-uniform coverage of drought characteristic space

## Sibling Projects

- `../SynHydro` — Kirsch-Nowak synthetic hydrology generator (Python package)
- `../NYCOptimization` — NYC reservoir optimization study (consumer of generated droughts)
- `../DiagnosticExperiment` or `../StochasticExploratoryExperiment` — ensemble simulation examples

## Rules

- Type hints and Google-style docstrings on all new code
- All experiments must be reproducible (config files, random seeds)
- Outputs are gitignored and regenerable
- Borg source files are licensed and gitignored
- **All `_pct` fields use 0-1 fractions. No exceptions.**
