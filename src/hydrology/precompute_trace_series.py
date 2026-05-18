"""Per-trace annual operational-outcome timeseries from Pywr-DRB output.

Used by Stage-09 magnitude-varying SA in the ``within_trace_percentile``
response form (Variant 2 of Hadjimichael et al. 2020 -- see
``magnitude_varying_sensitivity_analysis.py`` in
github.com/antonia-had/Magnitude_varying_sensitivity_analysis).

Hadjimichael's Variant 2 takes the within-SOW distribution of annual
shortage and feeds the τ-th percentile of *that distribution* to SALib
as the response variable Y_i(τ). The unit of analysis is the SOW; the
inner ensemble (10 stochastic realizations × 105 years) provides the
distribution that the percentile is taken over.

The MOEA-FIND analogue: each archive trace is a single 20-year
Pywr-DRB simulation. The analogue of the SOW's inner-distribution is
the trace's own annual sequence (20 values). For each percentile τ,
Y_i(τ) = np.percentile(annual_series_i, 100 * τ) -- a continuous
scalar per trace per percentile, no binary degeneracy, no synthetic
inner-ensemble fabrication. This module precomputes the per-trace ×
per-year matrix once so the SA driver can take percentiles cheaply.

Public surface:

- :func:`compute_annual_nyc_min_storage_frac` -- 3-reservoir aggregate
  storage / NYC_TOTAL_CAPACITY, calendar-year minimum per (year, trace).

Output schema: parquet with float index ``year_index`` (0..T-1) and
string columns = trace realization ids; values are the annual statistic.
The orientation matches the Stage-06 metric-bank index (string ids), so
downstream SA loaders can join on realization id directly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import h5py
import numpy as np
import pandas as pd

# Match src.pywrdrb.satisficing_metrics.NYC_STORAGE_CAPACITIES (same source file).
NYC_STORAGE_CAPACITIES = {
    "cannonsville": 95700.0,
    "pepacton": 140200.0,
    "neversink": 34900.0,
}
NYC_TOTAL_CAPACITY = float(sum(NYC_STORAGE_CAPACITIES.values()))


def _decode_time(time_dataset) -> pd.DatetimeIndex:
    """Decode the hdf5 ``time`` dataset (bytes) to a DatetimeIndex."""
    raw = time_dataset[:]
    if raw.dtype.kind == "S" or raw.dtype.kind == "O":
        as_str = [t.decode("utf-8") if isinstance(t, bytes) else str(t)
                  for t in raw]
    else:
        as_str = [str(t) for t in raw]
    return pd.to_datetime(as_str)


def compute_annual_nyc_min_storage_frac(
    pywrdrb_output_path: Path,
    *,
    realization_ids: Optional[list] = None,
) -> pd.DataFrame:
    """Annual minimum NYC reservoir storage fraction per trace.

    For each calendar year and each trace, the minimum daily aggregated
    NYC storage fraction (sum of cannonsville + pepacton + neversink
    storage, divided by their combined nameplate capacity ~270.8 BG).

    Args:
        pywrdrb_output_path: Path to ``simulations/pywrdrb_output.hdf5``
            written by Stage 06.
        realization_ids: Optional list of trace ids to use as column
            labels. Defaults to ``["0", "1", ..., "N-1"]``, matching the
            int->str convention used by the Stage-06 metric bank.

    Returns:
        DataFrame indexed by ``year`` (calendar year, int) with one column
        per realization id. Values are the annual minimum storage
        fraction in [0, 1].
    """
    pywrdrb_output_path = Path(pywrdrb_output_path)
    with h5py.File(pywrdrb_output_path, "r") as f:
        time_index = _decode_time(f["time"])
        n_traces = f["reservoir_cannonsville"].shape[1]
        # Sum daily storage across the three NYC reservoirs (units MG).
        # Reading the full (T, N) arrays at once; for our 7305 x 3308 grid
        # this is ~580 MB total which fits comfortably.
        agg = (
            f["reservoir_cannonsville"][:]
            + f["reservoir_pepacton"][:]
            + f["reservoir_neversink"][:]
        )
    frac = agg / NYC_TOTAL_CAPACITY  # daily fraction per (day, trace)

    if realization_ids is None:
        realization_ids = [str(i) for i in range(n_traces)]
    if len(realization_ids) != n_traces:
        raise ValueError(
            f"realization_ids has {len(realization_ids)} entries but "
            f"hdf5 has {n_traces} trace columns."
        )

    df_daily = pd.DataFrame(
        frac, index=time_index, columns=realization_ids,
    )
    df_annual = df_daily.groupby(df_daily.index.year).min()
    df_annual.index.name = "year"
    return df_annual
