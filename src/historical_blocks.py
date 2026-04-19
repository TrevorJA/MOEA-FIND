"""Multi-block resampling of a historical monthly record.

Synthetic traces produced by MOEA-FIND have a fixed length of ``T`` water
years. Comparing a single 1D historical record (typically much longer
than ``T``) against the synthetic ensemble biases the comparison: the
historical record has more months to sample its extreme tails than any
individual ``T``-year synthetic trace, so the historical FDC naturally
dips lower at the 99 %+ exceedance band.

The fair comparator is the *distribution* of ``T``-year blocks drawn
from the historical record. For ``N_hist`` historical water years and
``T`` target years, stride-``s`` resampling produces
``floor((N_hist - T) / s) + 1`` overlapping blocks. Each block is the
same length as a synthetic trace, so per-block FDC / ACF / seasonal-cycle
statistics are directly comparable to per-trace statistics.

The functions below deliberately do not assume pandas; callers that
need a DatetimeIndex should build it externally.
"""

from __future__ import annotations

from typing import List, Sequence

import numpy as np


def resample_historical_blocks(
    monthly_1d: np.ndarray,
    T_years: int,
    stride: int = 1,
) -> List[np.ndarray]:
    """Return overlapping T-year windows of a historical monthly record.

    Args:
        monthly_1d: 1D array of monthly flows, length ``N_hist * 12``.
            The array is assumed to already be water-year-aligned (i.e.
            its first month is October of some year).
        T_years: Block length in water years. Each block has length
            ``T_years * 12``.
        stride: Year stride between consecutive blocks. ``stride=1`` gives
            the maximum number of overlapping blocks; larger strides
            reduce correlation between blocks at the cost of fewer
            samples.

    Returns:
        List of 1D arrays, each of length ``T_years * 12``. Order is
        chronological (earliest block first).

    Raises:
        ValueError: if the record is shorter than ``T_years`` or if
            ``stride < 1``.
    """
    x = np.asarray(monthly_1d, dtype=float)
    if x.ndim != 1:
        raise ValueError(f"monthly_1d must be 1D, got shape {x.shape}")
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")

    n_months = len(x)
    n_hist_years = n_months // 12
    if n_hist_years < T_years:
        raise ValueError(
            f"historical record has {n_hist_years} water years, "
            f"need at least T_years={T_years} to form one block"
        )

    months_per_block = T_years * 12
    step = stride * 12
    blocks: List[np.ndarray] = []
    start = 0
    while start + months_per_block <= n_months:
        blocks.append(x[start : start + months_per_block].copy())
        start += step
    return blocks


def resample_historical_blocks_2d(
    monthly_2d: np.ndarray,
    T_years: int,
    stride: int = 1,
) -> List[np.ndarray]:
    """Same as :func:`resample_historical_blocks` but returns 2D blocks.

    Args:
        monthly_2d: Historical monthly flows reshaped to (n_hist_years, 12),
            water-year-aligned (column 0 is October).
        T_years: Block length in water years.
        stride: Year stride between consecutive blocks.

    Returns:
        List of 2D arrays of shape ``(T_years, 12)``.
    """
    x = np.asarray(monthly_2d, dtype=float)
    if x.ndim != 2 or x.shape[1] != 12:
        raise ValueError(
            f"monthly_2d must have shape (n_years, 12), got {x.shape}"
        )
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")

    n_hist_years = x.shape[0]
    if n_hist_years < T_years:
        raise ValueError(
            f"historical record has {n_hist_years} water years, "
            f"need at least T_years={T_years} to form one block"
        )

    blocks: List[np.ndarray] = []
    start = 0
    while start + T_years <= n_hist_years:
        blocks.append(x[start : start + T_years, :].copy())
        start += stride
    return blocks


def compute_historical_block_chars(
    monthly_1d: np.ndarray,
    T_years: int,
    ssi_calc,
    objective_keys: Sequence[str],
    stride: int = 1,
    start_date: str = "2100-01-01",
) -> np.ndarray:
    """Drought characteristics per historical T-year block.

    Extracts all overlapping T-year windows of ``monthly_1d`` via
    :func:`resample_historical_blocks`, transforms each block through
    ``ssi_calc`` (a SynHydro SSI calculator already fitted on the full
    historical record — **not** per-block, to mirror the DD-11 lock-in
    that drought classification uses the same calibrated distribution
    for all traces), and computes the drought characteristic vector for
    every block. The result is a ``(n_blocks, K)`` array in the same
    objective-key order Borg sees — ready to drop on top of the Pareto
    drought-space scatter as a historical cloud.

    Args:
        monthly_1d: 1D historical monthly flows, water-year-aligned.
        T_years: Block length (should match synthetic trace length).
        ssi_calc: Prefitted SynHydro SSI calculator (see
            :func:`src.objectives.make_ssi_calculator` +
            :func:`src.objectives.compute_ssi` with
            ``reference_flows=monthly_1d``).
        objective_keys: Drought characteristic keys to extract, same
            order as the MOEA objective vector.
        stride: Year stride between blocks (default 1 for maximum
            overlap).
        start_date: Dummy start date for the SSI's pandas Series
            (SSI is a stationary transform on fitted monthly
            distributions, so any well-formed date works).

    Returns:
        Array of shape ``(n_blocks, len(objective_keys))``.
    """
    from src.objectives import (
        flows_to_series,
        compute_ssi_drought_characteristics,
    )

    blocks = resample_historical_blocks(monthly_1d, T_years, stride)
    rows: List[list] = []
    for blk in blocks:
        ssi = ssi_calc.transform(flows_to_series(blk, start_date=start_date))
        chars = compute_ssi_drought_characteristics(ssi)
        rows.append([float(chars.get(k, 0.0)) for k in objective_keys])
    return np.asarray(rows, dtype=float)
