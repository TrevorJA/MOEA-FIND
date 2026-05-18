"""First-drought-event window extraction for the windowed pywrdrb pipeline.

MOEA-FIND optimises the characteristics of the **first critical SSI-3
drought** in each synthetic trace. Propagating the full 10-year trace
through Pywr-DRB then dilutes both the (expensive) simulation and the
downstream drought-performance metrics with non-drought years that the
optimisation never targeted. This module trims each scenario to just the
first-event window plus a justified buffer.

Window rule (2026-05-15 design decision)
----------------------------------------
For the first critical SSI-3 event with onset month ``o`` (event
``start``) and recovery month ``r`` (event ``end``)::

    window = [ o - LEAD_MONTHS , r + LAG_MONTHS ]   (clamped to the trace)

* ``LEAD_MONTHS = 36`` — the lead buffer serves a *dual* purpose:
  drought context **and** Pywr-DRB reservoir spin-up. DRB reservoirs
  (Cannonsville/Pepacton/Neversink + the STARFIT reservoirs) have long
  storage memory, so the windowed simulation must start far enough
  before onset for storage to equilibrate from the model's default
  initial state. 36 months (3 yr) is the chosen, defensible spin-up;
  it is asymmetric with the trailing buffer by design.
* ``LAG_MONTHS = 12`` — one year past recovery captures the
  post-drought refill / recovery dynamics that the performance metrics
  care about.

Edge handling: the window is clamped to ``[0, n_months]``. The realised
lead/lag are recorded per scenario; ``short_spinup`` flags scenarios
whose clamped lead is below ``LEAD_MONTHS`` (first event too early in
the trace) so downstream analysis / the manuscript can sensitivity-check
or exclude them rather than silently trusting an under-spun reservoir
state.

Resolution: bounds are computed at **monthly** resolution (SSI-3 is
monthly). :func:`window_dates` converts them to daily simulation
start/end dates for Pywr-DRB given the trace's calendar start.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Optional

import numpy as np
import pandas as pd

#: Lead buffer / reservoir spin-up, in months (see module docstring).
LEAD_MONTHS: int = 36
#: Trailing buffer past recovery, in months.
LAG_MONTHS: int = 12
#: Critical-drought SSI threshold (matches SynHydro / the stage-04 objective
#: definition: an event is "critical" once SSI <= -1).
CRITICAL_SSI: float = -1.0


@dataclass(frozen=True)
class DroughtWindow:
    """First-critical-event window for one scenario (monthly indices).

    Indices are 0-based month offsets into the trace. ``end_idx`` is
    *exclusive* (Python-slice convention): ``trace[start_idx:end_idx]``.
    """

    onset_idx: int          #: month index of first critical event start
    recovery_idx: int       #: month index of first critical event end
    start_idx: int          #: window start (clamped)  = max(0, onset-LEAD)
    end_idx: int            #: window end exclusive (clamped)
    realized_lead: int      #: onset_idx - start_idx (<= LEAD_MONTHS)
    realized_lag: int       #: (end_idx-1) - recovery_idx (<= LAG_MONTHS)
    n_window_months: int    #: end_idx - start_idx
    short_spinup: bool      #: realized_lead < LEAD_MONTHS (under-spun)

    def to_dict(self) -> Dict:
        return {k: (int(v) if isinstance(v, (np.integer,)) else v)
                for k, v in asdict(self).items()}


def compute_window(
    onset_idx: int,
    recovery_idx: int,
    n_months: int,
    *,
    lead_months: int = LEAD_MONTHS,
    lag_months: int = LAG_MONTHS,
) -> DroughtWindow:
    """Pure window arithmetic: clamp [onset-lead, recovery+lag] to the trace.

    Separated from SynHydro event location so it is unit-testable
    without SynHydro (the suite mocks ``synhydro`` at import time).
    """
    onset_idx = int(max(0, min(onset_idx, n_months - 1)))
    recovery_idx = int(max(onset_idx, min(recovery_idx, n_months - 1)))
    start_idx = max(0, onset_idx - int(lead_months))
    end_idx = min(n_months, recovery_idx + int(lag_months) + 1)
    realized_lead = onset_idx - start_idx
    realized_lag = (end_idx - 1) - recovery_idx
    return DroughtWindow(
        onset_idx=onset_idx,
        recovery_idx=recovery_idx,
        start_idx=start_idx,
        end_idx=end_idx,
        realized_lead=realized_lead,
        realized_lag=realized_lag,
        n_window_months=end_idx - start_idx,
        short_spinup=realized_lead < int(lead_months),
    )


def first_event_window(
    ssi_series: pd.Series,
    *,
    lead_months: int = LEAD_MONTHS,
    lag_months: int = LAG_MONTHS,
    end_drought_threshold_months: int = 3,
) -> Optional[DroughtWindow]:
    """Locate the first critical SSI-3 event and its buffered window.

    Args:
        ssi_series: SSI-3 ``pd.Series`` (monthly, DatetimeIndex) for the
            scenario — the same series used to compute the stage-04
            first-event objectives, so the window is consistent with
            what was optimised.
        lead_months / lag_months: buffer rule (see module docstring).
        end_drought_threshold_months: recovery hysteresis passed through
            to SynHydro's event extractor (production default 3).

    Returns:
        :class:`DroughtWindow`, or ``None`` if the trace contains no
        critical SSI-3 event (should not happen for stage-04-feasible
        scenarios, which pass the ">=1 critical event" hard constraint,
        but callers must handle ``None`` defensively).
    """
    from synhydro.droughts.ssi import get_drought_metrics

    dm = get_drought_metrics(
        ssi_series, end_drought_threshold_months=end_drought_threshold_months
    )
    if dm is None or len(dm) == 0:
        return None

    # First *critical* event in chronological order. get_drought_metrics
    # already records only critical droughts (severity reaches <= -1), so
    # the chronologically-first row is the first critical event.
    dm_sorted = dm.sort_values("start")
    first = dm_sorted.iloc[0]

    index = ssi_series.index
    n_months = len(index)

    def _month_pos(ts) -> int:
        ts = pd.Timestamp(ts)
        # Position of the event timestamp within the SSI monthly index.
        # searchsorted is robust to month-start vs month-any timestamps.
        pos = int(index.searchsorted(ts))
        return max(0, min(pos, n_months - 1))

    onset_idx = _month_pos(first["start"])
    recovery_idx = _month_pos(first["end"])
    return compute_window(
        onset_idx, recovery_idx, n_months,
        lead_months=lead_months, lag_months=lag_months,
    )


def window_dates(
    win: DroughtWindow,
    trace_start: str | pd.Timestamp,
    freq: str = "MS",
) -> Dict[str, pd.Timestamp]:
    """Convert monthly window indices to Pywr-DRB sim start/end dates.

    The simulation runs daily; ``sim_start`` is the first day of the
    window's first month and ``sim_end`` is the last day of the window's
    last month, so the daily Pywr-DRB run spans exactly the windowed
    months.

    Returns ``{"sim_start": Timestamp, "sim_end": Timestamp}``.
    """
    months = pd.date_range(start=pd.Timestamp(trace_start), periods=1, freq=freq)
    base = months[0]
    sim_start = (base + pd.DateOffset(months=int(win.start_idx))).normalize()
    # end_idx is exclusive in months; last included month = end_idx-1.
    last_month_start = base + pd.DateOffset(months=int(win.end_idx) - 1)
    sim_end = (last_month_start + pd.DateOffset(months=1)
               - pd.Timedelta(days=1)).normalize()
    return {"sim_start": sim_start, "sim_end": sim_end}
