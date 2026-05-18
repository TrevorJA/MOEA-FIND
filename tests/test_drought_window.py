"""Window-logic tests for src.metrics.drought_window.

The suite mocks ``synhydro`` at import time (see tests/conftest.py), so
these target the pure, SynHydro-free pieces: ``compute_window`` (the
clamp/buffer arithmetic) and ``window_dates`` (index→date conversion).
The thin SynHydro event-location wrapper ``first_event_window`` is
exercised in real SLURM runs, not here.
"""

import sys
from unittest.mock import MagicMock

import pandas as pd
import pytest

from src.metrics.drought_window import (
    DroughtWindow,
    compute_window,
    first_event_window,
    window_dates,
    LEAD_MONTHS,
    LAG_MONTHS,
)


def test_short_spinup_when_event_early_lead_clamped():
    w = compute_window(onset_idx=6, recovery_idx=11, n_months=180)
    assert w.start_idx == 0                 # max(0, 6-36)
    assert w.realized_lead == 6             # < LEAD_MONTHS
    assert w.short_spinup is True
    assert w.end_idx == 11 + LAG_MONTHS + 1
    assert w.n_window_months == w.end_idx - w.start_idx


def test_full_buffer_when_event_mid_trace():
    w = compute_window(onset_idx=80, recovery_idx=87, n_months=240)
    assert w.start_idx == 80 - LEAD_MONTHS
    assert w.realized_lead == LEAD_MONTHS
    assert w.short_spinup is False
    assert w.realized_lag == LAG_MONTHS


def test_trailing_clamp_near_trace_end():
    w = compute_window(onset_idx=100, recovery_idx=109, n_months=120)
    assert w.end_idx == 120                 # clamped to trace end
    assert w.realized_lag == (120 - 1) - 109
    assert w.realized_lead == LEAD_MONTHS


def test_recovery_before_onset_is_coerced():
    w = compute_window(onset_idx=50, recovery_idx=40, n_months=120)
    assert w.recovery_idx >= w.onset_idx    # coerced, no negative span
    assert w.n_window_months > 0


def test_window_dates_span_exact_months():
    w = DroughtWindow(
        onset_idx=40, recovery_idx=48, start_idx=4, end_idx=61,
        realized_lead=36, realized_lag=12, n_window_months=57,
        short_spinup=False,
    )
    d = window_dates(w, "1950-10-01")
    assert d["sim_start"] == pd.Timestamp("1951-02-01")  # month 4
    assert d["sim_end"].is_month_end
    assert d["sim_end"] > d["sim_start"]


@pytest.mark.skipif(
    isinstance(sys.modules.get("synhydro"), MagicMock),
    reason="synhydro is mocked in the test suite; first_event_window "
           "needs real get_drought_metrics (covered by SLURM runs).",
)
def test_first_event_window_integration():
    import numpy as np
    n = 240
    vals = np.full(n, 0.2)
    vals[80:88] = -1.6
    idx = pd.date_range("1950-10-01", periods=n, freq="MS")
    w = first_event_window(pd.Series(vals, index=idx))
    assert w is not None and w.short_spinup is False
