"""Safety-net tests for src.pywrdrb.satisficing_metrics pure helpers.

Phase 0 regression net for the src/ -> src/pywrdrb/ move. Exercises the
pure deficit/storage/level extractors (no Pywr-DRB HDF5 needed).
"""

import numpy as np
import pandas as pd

from src.pywrdrb.satisficing_metrics import (
    _longest_run,
    _hashimoto,
    _ffmp_exposure,
    _storage_minima,
    _safe_first_column,
    NYC_STORAGE_CAPACITIES,
    NYC_TOTAL_CAPACITY,
)


class TestLongestRun:
    def test_all_false(self):
        assert _longest_run(np.zeros(5, dtype=bool)) == 0

    def test_all_true(self):
        assert _longest_run(np.ones(5, dtype=bool)) == 5

    def test_split_runs(self):
        m = np.array([1, 1, 0, 1, 1, 1, 0], dtype=bool)
        assert _longest_run(m) == 3


class TestSafeFirstColumn:
    def test_picks_first_present(self):
        df = pd.DataFrame({"b": [1], "a": [2]})
        s = _safe_first_column(df, ("a", "b"))
        assert s is not None and s.iloc[0] == 2

    def test_none_when_absent(self):
        df = pd.DataFrame({"x": [1]})
        assert _safe_first_column(df, ("a", "b")) is None


class TestHashimoto:
    def test_none_inputs_return_nan_dict(self):
        out = _hashimoto(None, None)
        assert set(out) == {
            "reliability", "resiliency", "vulnerability",
            "max_event_duration_days", "event_count",
        }
        assert np.isnan(out["reliability"])

    def test_perfect_reliability(self):
        flow = pd.Series(np.full(50, 100.0))
        target = pd.Series(np.full(50, 100.0))
        out = _hashimoto(flow, target)
        assert out["reliability"] == 1.0
        assert out["vulnerability"] == 0.0
        assert out["event_count"] == 0.0

    def test_metrics_bounded_and_monotone_in_deficit(self):
        target = pd.Series(np.full(100, 100.0))
        mild = pd.Series(np.where(np.arange(100) < 10, 90.0, 100.0))
        severe = pd.Series(np.where(np.arange(100) < 40, 50.0, 100.0))
        m_mild = _hashimoto(mild, target)
        m_sev = _hashimoto(severe, target)
        for m in (m_mild, m_sev):
            assert 0.0 <= m["reliability"] <= 1.0
        # More/deeper deficit -> lower reliability, larger vulnerability.
        assert m_sev["reliability"] < m_mild["reliability"]
        assert m_sev["vulnerability"] >= m_mild["vulnerability"]


class TestFFMPExposure:
    def test_none_returns_all_nan(self):
        out = _ffmp_exposure(None)
        assert np.isnan(out["max_level"])

    def test_level_counts(self):
        lv = pd.DataFrame({"nyc": [0, 1, 3, 6, 6, 2]})
        out = _ffmp_exposure(lv)
        assert out["max_level"] == 6
        assert out["days_drought_any"] == 5
        assert out["days_L6"] == 2
        assert out["first_L6_day"] == 3


class TestStorageMinima:
    def test_none_returns_nan_keys(self):
        out = _storage_minima(None)
        assert np.isnan(out["nyc_min_storage_frac"])

    def test_fraction_in_unit_range_and_nan_free(self):
        n = 30
        df = pd.DataFrame({
            r: np.linspace(cap, cap * 0.4, n)
            for r, cap in NYC_STORAGE_CAPACITIES.items()
        })
        out = _storage_minima(df)
        frac = out["nyc_min_storage_frac"]
        assert 0.0 <= frac <= 1.0
        assert not np.isnan(frac)
        # Min at the end of the linspace: ~0.4 of total capacity.
        assert frac == np.float64(
            sum(NYC_STORAGE_CAPACITIES.values()) * 0.4 / NYC_TOTAL_CAPACITY
        )
