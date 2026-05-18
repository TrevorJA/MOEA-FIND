"""Shared SSI Tier-A event-summary extractor.

Both :func:`src.metrics.objectives.compute_ssi_drought_characteristics`
and :func:`src.metrics.extended.compute_ssi_event_metrics` derived the
same 10 Tier-A event statistics from a ``get_drought_metrics`` table.
That logic now lives here once; the two callers delegate to it (the
former for its core block, then layering first-event/peak extras on
top; the latter as a thin suffixed wrapper).

The extraction is verbatim — see ``tests/golden/ssi_common_gate.py``
for the behavior-preservation gate run against real SynHydro.
"""

from typing import Dict, Optional

import pandas as pd
from synhydro.droughts.ssi import get_drought_metrics

from src.metrics.objectives import TIME_IN_DROUGHT_THRESHOLD


def compute_ssi_tier_a(
    ssi_series: pd.Series,
    *,
    dm: Optional[pd.DataFrame] = None,
    suffix: str = "",
    end_drought_threshold_months: int = 3,
) -> Dict[str, float]:
    """The 10 standard SSI Tier-A event metrics from an SSI series.

    Args:
        ssi_series: SSI :class:`pandas.Series` (DatetimeIndex required).
        dm: Optional pre-computed ``get_drought_metrics`` table. When a
            caller already has it (objectives needs it for the
            first-event/peak extras), pass it here to avoid recomputing
            — the result is identical either way.
        suffix: Appended to every output key (``""`` reproduces the
            production SSI-3 names; ``"_ssi12"`` for the long-timescale
            candidate set).
        end_drought_threshold_months: Recovery hysteresis (months of
            consecutive SSI > 0 to terminate a critical drought).

    Returns:
        Dict with keys (each with ``suffix`` appended): ``frequency``,
        ``mean_duration``, ``max_duration``, ``mean_magnitude``,
        ``max_magnitude``, ``mean_severity``, ``worst_severity``,
        ``mean_avg_severity``, ``time_in_drought_fraction``, ``n_events``.
    """
    if dm is None:
        dm = get_drought_metrics(
            ssi_series,
            end_drought_threshold_months=end_drought_threshold_months,
        )

    valid = ssi_series.dropna()
    n_valid = len(valid)
    n_years = n_valid / 12.0 if n_valid > 0 else 0.0

    if n_valid > 0:
        time_in_drought = float(
            (valid <= TIME_IN_DROUGHT_THRESHOLD).sum() / n_valid
        )
    else:
        time_in_drought = 0.0

    if len(dm) == 0:
        return {
            f"frequency{suffix}": 0.0,
            f"mean_duration{suffix}": 0.0,
            f"max_duration{suffix}": 0.0,
            f"mean_magnitude{suffix}": 0.0,
            f"max_magnitude{suffix}": 0.0,
            f"mean_severity{suffix}": 0.0,
            f"worst_severity{suffix}": 0.0,
            f"mean_avg_severity{suffix}": 0.0,
            f"time_in_drought_fraction{suffix}": time_in_drought,
            f"n_events{suffix}": 0,
        }

    return {
        f"frequency{suffix}": float(len(dm) / n_years * 10) if n_years > 0 else 0.0,
        f"mean_duration{suffix}": float(dm["duration"].mean()),
        f"max_duration{suffix}": float(dm["duration"].max()),
        f"mean_magnitude{suffix}": float(dm["magnitude"].abs().mean()),
        f"max_magnitude{suffix}": float(dm["magnitude"].abs().max()),
        f"mean_severity{suffix}": float(dm["severity"].abs().mean()),
        f"worst_severity{suffix}": float(dm["severity"].abs().max()),
        f"mean_avg_severity{suffix}": float(dm["avg_severity"].abs().mean()),
        f"time_in_drought_fraction{suffix}": time_in_drought,
        f"n_events{suffix}": int(len(dm)),
    }
