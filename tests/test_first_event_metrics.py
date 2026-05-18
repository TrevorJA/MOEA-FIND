"""Tests for the first-event SSI-3 drought objective family (T=10y pivot).

Drives compute_ssi_drought_characteristics by monkey-patching the imported
get_drought_metrics symbol with a stub that returns controlled DataFrames,
so the tests run without SynHydro installed (conftest stubs synhydro).
"""

import numpy as np
import pandas as pd
import pytest

from src.metrics import objectives
from src.metrics.drought_metrics import (
    REGISTRY,
    PRESETS,
    AntiIdealRule,
    compute_anti_ideal,
    resolve_metric_set,
)


def _fake_ssi_series(n_months: int = 36) -> pd.Series:
    """Plausible SSI series with monthly DatetimeIndex."""
    idx = pd.date_range("2000-01-01", periods=n_months, freq="MS")
    return pd.Series(np.zeros(n_months), index=idx)


def _patched_get_drought_metrics(monkeypatch, dm_df: pd.DataFrame):
    """Replace src.metrics.objectives.get_drought_metrics with a stub returning dm_df."""
    monkeypatch.setattr(
        objectives, "get_drought_metrics", lambda *a, **kw: dm_df.copy()
    )


# ---------------------------------------------------------------------------
# Single-event extraction
# ---------------------------------------------------------------------------


class TestSingleEvent:
    def test_first_event_metrics_match_hand_computed(self, monkeypatch):
        # One event: starts 2001-04, peak 2001-07, ends 2001-10. Duration 7 mo.
        # severity (min SSI) = -2.3 → |severity| = 2.3
        # magnitude (sum SSI) = -8.5 → |magnitude| = 8.5
        dm = pd.DataFrame([{
            "start": pd.Timestamp("2001-04-01"),
            "end": pd.Timestamp("2001-10-01"),
            "duration": 7,
            "magnitude": -8.5,
            "severity": -2.3,
            "avg_severity": -1.214,
            "max_severity_date": pd.Timestamp("2001-07-01"),
            "recovery_period": 3,
        }])
        _patched_get_drought_metrics(monkeypatch, dm)

        chars = objectives.compute_ssi_drought_characteristics(_fake_ssi_series())

        assert chars["first_event_present"] == 1
        assert chars["first_event_duration"] == pytest.approx(7.0)
        assert chars["first_event_severity"] == pytest.approx(2.3)
        assert chars["first_event_magnitude"] == pytest.approx(8.5)
        assert chars["first_event_start_month"] == pytest.approx(4.0)
        assert chars["first_event_peak_month"] == pytest.approx(7.0)


# ---------------------------------------------------------------------------
# Zero-event behavior (empty dm) — drives the at-least-one-event constraint
# ---------------------------------------------------------------------------


class TestZeroEvent:
    def test_empty_dm_zero_first_event_keys(self, monkeypatch):
        empty = pd.DataFrame()
        _patched_get_drought_metrics(monkeypatch, empty)
        chars = objectives.compute_ssi_drought_characteristics(_fake_ssi_series())

        assert chars["first_event_present"] == 0
        assert chars["first_event_duration"] == 0.0
        assert chars["first_event_severity"] == 0.0
        assert chars["first_event_magnitude"] == 0.0
        assert chars["first_event_start_month"] == 0.0
        assert chars["first_event_peak_month"] == 0.0
        assert chars["n_events"] == 0


# ---------------------------------------------------------------------------
# Multi-event: pick the chronologically-first row
# ---------------------------------------------------------------------------


class TestMultiEvent:
    def test_picks_chronologically_first_event(self, monkeypatch):
        # Two events: build the dm in REVERSE chronological order to confirm
        # the sort_values("start") in the implementation actually runs.
        dm = pd.DataFrame([
            {
                "start": pd.Timestamp("2005-09-01"),
                "end": pd.Timestamp("2006-02-01"),
                "duration": 6,
                "magnitude": -12.0,
                "severity": -3.1,
                "avg_severity": -2.0,
                "max_severity_date": pd.Timestamp("2005-12-01"),
                "recovery_period": 3,
            },
            {
                "start": pd.Timestamp("2002-03-01"),
                "end": pd.Timestamp("2002-08-01"),
                "duration": 6,
                "magnitude": -7.2,
                "severity": -1.8,
                "avg_severity": -1.2,
                "max_severity_date": pd.Timestamp("2002-05-01"),
                "recovery_period": 3,
            },
        ])
        _patched_get_drought_metrics(monkeypatch, dm)
        chars = objectives.compute_ssi_drought_characteristics(_fake_ssi_series())

        # Should pick the 2002 event, not the 2005 event.
        assert chars["first_event_severity"] == pytest.approx(1.8)
        assert chars["first_event_magnitude"] == pytest.approx(7.2)
        assert chars["first_event_start_month"] == pytest.approx(3.0)
        assert chars["first_event_peak_month"] == pytest.approx(5.0)
        # Aggregate keys still cover all events.
        assert chars["n_events"] == 2
        assert chars["worst_severity"] == pytest.approx(3.1)
        assert chars["max_magnitude"] == pytest.approx(12.0)


# ---------------------------------------------------------------------------
# December/January wraparound — the start_month and peak_month are calendar
# months 1–12, so a December event reports 12 (not 0 or 13).
# ---------------------------------------------------------------------------


class TestCyclicMonthEdgeCases:
    def test_december_event(self, monkeypatch):
        dm = pd.DataFrame([{
            "start": pd.Timestamp("2003-12-01"),
            "end": pd.Timestamp("2004-04-01"),
            "duration": 5,
            "magnitude": -6.0,
            "severity": -1.5,
            "avg_severity": -1.2,
            "max_severity_date": pd.Timestamp("2004-01-01"),
            "recovery_period": 3,
        }])
        _patched_get_drought_metrics(monkeypatch, dm)
        chars = objectives.compute_ssi_drought_characteristics(_fake_ssi_series())

        assert chars["first_event_start_month"] == pytest.approx(12.0)
        assert chars["first_event_peak_month"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Registry & preset wiring
# ---------------------------------------------------------------------------


class TestRegistryAndPreset:
    def test_preset_resolves_to_five_metrics(self):
        preset = PRESETS["first_event_ssi3_t10"]
        assert len(preset) == 5
        metric_set = resolve_metric_set("first_event_ssi3_t10")
        assert len(metric_set) == 5
        names = [m.name for m in metric_set]
        assert names == [
            "first_event_duration",
            "first_event_severity",
            "first_event_magnitude",
            "first_event_start_month",
            "first_event_peak_month",
        ]

    def test_max_partners_anchor_to_all_events_keys(self):
        # Anti-ideal anchors first-event metrics against all-events historical
        # maxes (worst_severity, max_magnitude, max_duration), so re-running
        # against an isolated single-event historical dm produces D* values
        # equal to headroom × those all-events maxes.
        assert REGISTRY["first_event_severity"].max_partner == "worst_severity"
        assert REGISTRY["first_event_magnitude"].max_partner == "max_magnitude"
        assert REGISTRY["first_event_duration"].max_partner == "max_duration"

    def test_month_metrics_are_cyclic(self):
        for name in ("first_event_start_month", "first_event_peak_month"):
            m = REGISTRY[name]
            assert m.is_cyclic is True
            assert m.anti_ideal_rule is AntiIdealRule.CYCLIC_HEADROOM

    def test_continuous_metrics_use_headroom_times_max(self):
        for name in ("first_event_duration", "first_event_severity",
                     "first_event_magnitude"):
            m = REGISTRY[name]
            assert m.is_cyclic is False
            assert m.anti_ideal_rule is AntiIdealRule.HEADROOM_TIMES_MAX


# ---------------------------------------------------------------------------
# Anti-ideal: the all-events historical envelope drives D*
# ---------------------------------------------------------------------------


class TestAntiIdealAnchoring:
    def test_anti_ideal_uses_partner_keys(self):
        metric_set = resolve_metric_set("first_event_ssi3_t10")
        # hist_chars mimics what compute_ssi_drought_characteristics returns
        # on the FULL historical record (not just first event).
        hist_chars = {
            "max_duration": 18.0,
            "worst_severity": 2.7,
            "max_magnitude": 30.0,
        }
        D = compute_anti_ideal(metric_set, hist_chars, headroom=1.5)
        # Order matches PRESETS["first_event_ssi3_t10"]: duration, severity,
        # magnitude, start_month, peak_month.
        assert D[0] == pytest.approx(1.5 * 18.0)   # duration
        assert D[1] == pytest.approx(1.5 * 2.7)    # severity (peak depth)
        assert D[2] == pytest.approx(1.5 * 30.0)   # magnitude (cumulative)
        assert D[3] == pytest.approx(12.0 * 1.5)   # cyclic month
        assert D[4] == pytest.approx(12.0 * 1.5)   # cyclic month
