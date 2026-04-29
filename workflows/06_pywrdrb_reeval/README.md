# Stage 06 -- Pywr-DRB policy re-evaluation

## Purpose

Re-evaluate each Pareto drought scenario from stage 04 through Pywr-DRB
under the FFMP NYC operating policy and `constant_max` demand. Produces
a per-realization metric bank (FFMP exposure, NYC storage, flow targets,
delivery reliability) that feeds the satisficing / scenario-discovery
analysis in stage 07.

`policy_reeval.py` invokes `replay_pareto_to_multisite_monthly` from
`src.pywrdrb_bridge`, which preserves spatial correlation across the
multi-site ensemble by construction: every Pareto DV vector is replayed
through *shared* monthly Kirsch indexes, so all sites see the same
monthly resampling pattern. The legacy "stage 05" multi-site driver has
been folded into this step and is no longer needed as a separate stage.

## Drivers

| Driver | Purpose | Manuscript |
|---|---|---|
| [verify_drought_coverage.py](verify_drought_coverage.py) | Pre-flight pass/fail report on the stage 04 archive (Criteria 1-5).         | SI-I |
| [policy_reeval.py](policy_reeval.py)                     | Replay Pareto DVs through Pywr-DRB; metric bank + FFMP satisficing table.   | §7.3 main |

`verify_drought_coverage` should be run **before** committing a long
`policy_reeval` job: it confirms the input Pareto archive covers the
drought-characteristic space without gaps. `policy_reeval` is the
expensive step (hours on a full archive) and feeds stage 07.

## Compute / plot split

Compute drivers write only numerical artifacts under
`outputs/06_pywrdrb_reeval/<driver>/<src_slug>/`. Figures are produced
by the paired plotting drivers under [plots/](plots/) which read those
artifacts and write to
`figures/06_pywrdrb_reeval/<driver>/<src_slug>/`. Re-rendering a figure
never requires re-running the simulation.

`<src_slug>` is the parent-directory name of the input
`results.json` (e.g. `residual_T20_nfe200000_s42_constrained_cmdv_uniform_stad`),
so stage 06 outputs pair 1:1 with their stage 04 archive on disk.

## Slurm

| Compute slurm | Plotting slurm |
|---|---|
| `slurm/policy_reeval.slurm`           | -- (figures emitted by stage 07)                |
| `slurm/verify_drought_coverage.slurm` | `slurm/plots/verify_drought_coverage.slurm`     |

`policy_reeval` requests 3 nodes x 30 ntasks-per-node = 90 ranks (Stages
2-3 are MPI-parallel; Stages 1 and 4 run serially on rank 0).
`verify_drought_coverage` is 1 node, 4 cores.

All slurm scripts are self-contained -- the `PARETO_RESULTS` path is
baked at the top of each file. To re-evaluate or verify a different
stage 04 archive, edit the slurm script; do not pass `--export` arguments
via `sbatch`.

## Outputs

```
outputs/06_pywrdrb_reeval/
  policy_reeval/<src_slug>/
    config.json
    pywrdrb_inputs/{gage_flow_mgd.hdf5, catchment_inflow_mgd.hdf5,
                    predicted_inflows_mgd.hdf5}
    simulations/pywrdrb_output.hdf5
    results/{metric_bank.parquet, satisficing_table.parquet,
             drought_levels.npz}
  verify_drought_coverage/<src_slug>/
    config.json
    verification_report.json
    subsets.npz
    annual_means.npz
    hist_block_chars.npz
```

## Figures

```
figures/06_pywrdrb_reeval/
  verify_drought_coverage/<src_slug>/
    c1_drought_coverage.pdf
    c2_drought_subset_fdc.pdf
    c3_nominal_subset.pdf
    c4_annual_spread.pdf
    c5_drought_onset.pdf
```

`policy_reeval` emits no figures of its own -- the main satisficing
map (Fig 9) is rendered by stage 07's
`scenario_discovery_plots.py`, which consumes
`results/metric_bank.parquet` from this stage. Final manuscript PDFs
are regenerated from these working figures by
[workflows/99_manuscript_figures/make_figures.py](../99_manuscript_figures/make_figures.py).
