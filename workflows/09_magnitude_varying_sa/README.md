# Stage 09 -- Magnitude-Varying Sensitivity Analysis (MV-SA)

Adaptation of Hadjimichael et al. (2020) MV-SA to the MOEA-FIND
diagnostic question: at each percentile of an operational hazard
outcome (the *magnitude axis*, e.g. NYC minimum reservoir storage),
which **drought-hazard characteristics** drive the system response?

This is the methodological complement to Stage 07 scenario discovery:
where Stage 07 asks whether failure is *separable* in characteristic
space, Stage 09 asks which characteristic dimensions *dominate* at
each severity of operational stress.

## Method

The factor space is the optimized MOEA-FIND objective axes (drought
characteristics from the upstream archive, identical to Stage 08).
The magnitude axis is a single column from the Stage-06 metric bank.
For each percentile τ ∈ {0.05, 0.10, ..., 0.95}:

- **exceedance form** (default): build `I(τ) = 1{M ≤ M_quantile(τ)}`
  and run the SA method on `(X → I(τ))`. This measures which factors
  drive the *frequency* with which the magnitude axis crosses each
  percentile threshold (the direct analogue of Hadjimichael's
  "frequency response" at each shortage percentile).

- **conditional form**: restrict the sample to a window of width
  `window_frac × N` realizations centred on rank-τ, then run SA of
  the factors on a *secondary* outcome within that subset. Diagnostic
  for "given the system sits at this severity of the primary hazard,
  which factors drive the secondary outcome?"

A uniform-random *control factor* is appended (Hadjimichael et al.
2020) so the magnitude-varying noise floor is empirically visible at
every percentile.

## Files

- `run_mv_sa.py` -- compute driver. **MPI-parallel across percentiles**
  via `mpi4py`; works trivially serial when launched on one rank.
  Writes parquet + JSON under
  `outputs/09_magnitude_varying_sa/run_mv_sa/<slug>/`.
- `configs/delta_only.yaml` -- production preset (Delta only,
  19 percentiles, 50 bootstrap). Other methods are kept available in
  the engine but disabled here while the methodology is finalised.
- `configs/all_methods.yaml` -- multi-method preset (Delta + PAWN +
  RBD-FAST). Currently unused; available for later cross-method
  triangulation.
- `configs/delta.yaml` -- cheap dev preset (Delta only, coarse grid).
- `plots/run_mv_sa.py` -- plotting driver. Emits PDFs under
  `figures/09_magnitude_varying_sa/run_mv_sa/<slug>/`.
- `slurm/run_mv_sa.slurm` -- 19-rank MPI compute job (1 hr / 32 GB).
- `slurm/plots/run_mv_sa.slurm` -- 2-CPU/4-GB/30-minute plots job.

## Usage

```bash
# Production (19-rank MPI):
sbatch workflows/09_magnitude_varying_sa/slurm/run_mv_sa.slurm
sbatch workflows/09_magnitude_varying_sa/slurm/plots/run_mv_sa.slurm

# Dev / smoke (login node OK; trivially serial):
python workflows/09_magnitude_varying_sa/run_mv_sa.py \
    --bank  outputs/06_pywrdrb_reeval/policy_reeval/<src>/results/metric_bank.csv \
    --chars outputs/04_moea_find_single_site/run_moea_find/<src>/results.json \
    --config workflows/09_magnitude_varying_sa/configs/delta.yaml

# Multi-rank MPI from CLI:
mpirun -np 19 python workflows/09_magnitude_varying_sa/run_mv_sa.py \
    --bank  ... --chars ... \
    --config workflows/09_magnitude_varying_sa/configs/delta_only.yaml
```

## Output schema

`mv_sa_<method>.parquet` is long-form with columns:

| column | meaning |
| --- | --- |
| `percentile` | τ value in (0, 1) |
| `method` | `delta` / `pawn` / `rbd_fast` |
| `factor` | factor name (D characteristic or `control_uniform`) |
| `headline_index` | base method index at this τ |
| `ci_lo`, `ci_hi` | bootstrap CI on the headline index |
| `full_rank` | factor rank at this τ (1 = most sensitive) |
| `median_rank` | bootstrap median rank |
| `rank_iqr_lo`, `rank_iqr_hi` | bootstrap IQR of factor rank |
| `rank_spearman_median` | per-τ rank-Spearman vs full sample |
| `n_used` | rows entering this slice (full N or window subset) |
| `threshold` | M-quantile value at this τ (for axis annotation) |
| `axis_column` | metric-bank column used as M |
