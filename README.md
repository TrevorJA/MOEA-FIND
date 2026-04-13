# MOEA-FIND: Multi-Objective Search-Based Synthetic Drought Generation

Couples Borg MOEA with block bootstrap streamflow generation to produce ensembles of synthetic traces with structured coverage of drought characteristic space.

## Method

Existing search-based generators (Borgomeo et al. 2015; Zaniolo et al. 2023) use single-objective optimization to match one target drought at a time. MOEA-FIND formulates generation as a multi-objective problem:

- **Objectives:** Drought frequency, intensity, duration (relative to historical), plus an L1 Manhattan norm from the anti-ideal point
- **Decision variables:** Bootstrap sampling indices for the Kirsch-Nowak generator
- **Result:** The Pareto front is an ensemble of synthetic traces with near-uniform coverage of drought characteristic space

The Manhattan norm objective, combined with Borg's epsilon-dominance archiving, produces approximately uniform spacing along the Pareto hyperplane when projected to the drought characteristic dimensions.

## Installation

```bash
pip install -e .
```

Requires: `numpy`, `scipy`, `matplotlib`, `synhydro` (sibling package)

Borg MOEA binaries are not included (licensed). Place `borg.py` and compiled binaries in `lib/`.

## Usage

```bash
# Analytic proof of concept (2D/3D)
python experiments/proof_of_concept/run_analytic_poc.py

# Kirsch bootstrap with SSI objectives
python experiments/proof_of_concept/run_kirsch_poc.py --nfe 5000 --compare

# Library-and-subsample baseline
python experiments/kirsch_ensemble/run_library_baseline.py --n-traces 200
```
