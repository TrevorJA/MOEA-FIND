"""Experiment data loading & historical drought characterization.

Split out of the former ``src.experiment_utils`` god module: this file
owns USGS data prep and the historical-block characterization helpers
used to anchor anti-ideal points. No MOEA orchestration lives here.
"""

import json
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from src.io_paths.data import load_usgs_daily, daily_to_monthly
from src.metrics.objectives import (
    compute_ssi,
    compute_ssi_drought_characteristics,
)


def prepare_data(cache_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load historical USGS data and return (monthly_2d, monthly_1d).

    Args:
        cache_dir: Directory for caching downloaded USGS data.

    Returns:
        Tuple of (monthly_2d, monthly_1d) where:
        - monthly_2d: shape (n_years, 12)
        - monthly_1d: flattened historical flows
    """
    daily = load_usgs_daily(cache_dir=cache_dir)
    monthly = daily_to_monthly(daily)

    # Align to water year (October-September)
    first_oct = monthly.index[monthly.index.month == 10][0]
    last_sep = monthly.index[monthly.index.month == 9][-1]
    monthly = monthly[first_oct:last_sep]

    n_years = len(monthly) // 12
    monthly_values = monthly.values[:n_years * 12]
    monthly_2d = monthly_values.reshape(n_years, 12)

    print(f"Historical: {n_years} water years, mean={monthly_2d.mean():.1f} cfs")
    return monthly_2d, monthly_values


def compute_historical_ssi_chars(
    monthly_1d,
    timescale: int,
    **kwargs,
) -> Tuple:
    """Compute SSI and drought characteristics for historical data.

    Args:
        monthly_1d: 1D array or Series of monthly flows.
        timescale: SSI accumulation period (1, 3, 6, 12).

    Returns:
        Tuple of (ssi_series, ssi_calculator, drought_characteristics_dict).
    """
    ssi, calc = compute_ssi(monthly_1d, timescale=timescale)
    chars = compute_ssi_drought_characteristics(ssi)
    return ssi, calc, chars


def extract_pareto_maxes(
    results_json_path: Path,
    objective_keys,
) -> Dict[str, float]:
    """Read a prior ``results.json`` and return ``{name: max_D_j}`` per objective.

    Used to drive a Pareto-based anti-ideal placement on subsequent runs
    via :func:`src.experiment.anti_ideal.compute_ssi_anti_ideal`'s
    ``feasible_maxes`` argument.

    Args:
        results_json_path: Path to a prior ``results.json``.
        objective_keys: Either a tuple of metric names or a tuple of
            :class:`src.metrics.drought_metrics.DroughtMetric` instances.

    Raises:
        ValueError: ``objective_keys`` don't match the reference file's
            objective ordering. Mismatched objective sets would silently
            mis-scale ``D*``.
    """
    from src.metrics.drought_metrics import metric_names, resolve_metric_set

    metric_set = resolve_metric_set(objective_keys)
    names = metric_names(metric_set)

    payload = json.loads(Path(results_json_path).read_text())
    ref_keys = tuple(payload.get("objective_keys", ()))
    if names != ref_keys:
        raise ValueError(
            f"objective_keys mismatch: requested {names} but "
            f"{results_json_path} was run with {ref_keys}. Rebuild the "
            f"reference run with matching objectives, or drop the --anti-ideal-"
            f"reference flag to fall back to historical max."
        )
    dm = np.asarray(payload.get("drought_metrics", []), dtype=float)
    if dm.size == 0:
        raise ValueError(
            f"{results_json_path} has no drought_metrics (empty Pareto); "
            f"cannot derive feasible maxes."
        )
    return {n: float(dm[:, j].max()) for j, n in enumerate(names)}


# ---------------------------------------------------------------------------
# Short-block evaluation helpers (T=1 raw-flow path, DD-15)
# ---------------------------------------------------------------------------


def compute_historical_short_block_chars(
    monthly_2d: np.ndarray,
    monthly_1d: np.ndarray,
) -> dict:
    """Compute short-block metrics on each historical water year; return per-metric max.

    Used to set anti-ideal coordinates for raw-flow objectives at T=1.
    Each row of ``monthly_2d`` is one water year (12 months, Oct–Sep ordering).
    """
    from src.metrics.extended import FullRecordRefs
    from src.metrics.short_block import compute_short_block_metrics, SHORT_BLOCK_METRIC_NAMES

    refs = FullRecordRefs.from_full_record(monthly_1d)
    all_chars = []
    for yr_flows in monthly_2d:
        c = compute_short_block_metrics(yr_flows, ssi3_full_window_series=None, refs=refs)
        all_chars.append(c)
    maxes: dict = {}
    for key in SHORT_BLOCK_METRIC_NAMES:
        vals = [c.get(key, 0.0) for c in all_chars]
        maxes[key] = float(np.max(vals)) if vals else 0.0
    return maxes


def make_short_block_chars_fn(monthly_1d: np.ndarray):
    """Return an evaluate-compatible chars function for short-block metrics.

    Returns a flat dict combining:

    * Legacy 17-metric ``compute_short_block_metrics`` output
      (``djf_total_neg``, ``summer_recession``, ``aug_zscore``,
      ``ond_total_neg``, etc.) — preserved for backward-compatibility
      with archives created before the DD-15c reformulation.
    * 48-metric ``compute_candidate_bounded_metrics`` output (DD-15c
      Phase 1) — sub-annual log-space window summaries × {Mapping G,
      Mapping E}, all bounded in ``[0, 1)`` with ``D* = 1.0``. The
      production preset ``short_block_drb_v2`` selects K=4 keys from
      this set.

    Both reference structures (``FullRecordRefs`` and
    ``BoundedFamilyRefs``) are fitted once on the full historical
    record at factory construction time and reused on every evaluation.
    """
    from src.metrics.extended import (
        BoundedFamilyRefs,
        FullRecordRefs,
    )
    from src.metrics.short_block import (
        compute_candidate_bounded_metrics,
        compute_short_block_metrics,
    )

    refs = FullRecordRefs.from_full_record(monthly_1d)
    bounded_refs = BoundedFamilyRefs.from_full_record(monthly_1d)

    def _fn(synthetic_1d: np.ndarray, synthetic_2d: np.ndarray) -> dict:
        legacy = compute_short_block_metrics(
            synthetic_1d, ssi3_full_window_series=None, refs=refs,
        )
        bounded = compute_candidate_bounded_metrics(synthetic_1d, bounded_refs)
        legacy.update(bounded)
        return legacy

    return _fn
