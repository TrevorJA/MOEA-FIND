"""Unit tests for src.constraints — pure numpy, no SynHydro required."""

import numpy as np
import pytest

from src.constraints import (
    annual_stats_constraint,
    autocorrelation_constraint,
    seasonal_cycle_constraint,
)


# ---------------------------------------------------------------------------
# autocorrelation_constraint
# ---------------------------------------------------------------------------

class TestAutocorrelationConstraint:
    def test_identical_series_no_violation(self):
        rng = np.random.default_rng(0)
        x = rng.standard_normal(120)
        assert autocorrelation_constraint(x, x, tolerance=0.3) == pytest.approx(0.0)

    def test_positive_violation_for_constant_vs_random(self):
        rng = np.random.default_rng(1)
        random_series = rng.standard_normal(120)
        # Constant series has autocorr of NaN/0; random series has nonzero lag-1 autocorr
        # A very smooth (highly autocorrelated) series vs a random one should exceed tolerance
        smooth = np.sin(np.linspace(0, 4 * np.pi, 120))  # lag-1 autocorr ≈ 0.99
        violation = autocorrelation_constraint(smooth, random_series, tolerance=0.3)
        assert violation > 0.0

    def test_tight_tolerance_raises_violation(self):
        rng = np.random.default_rng(2)
        x = rng.standard_normal(120)
        y = rng.standard_normal(120)
        # With tolerance=0, any non-identical autocorrelation produces a violation
        v0 = autocorrelation_constraint(x, y, tolerance=0.0)
        v1 = autocorrelation_constraint(x, y, tolerance=0.3)
        assert v0 >= v1

    def test_violation_is_non_negative(self):
        rng = np.random.default_rng(3)
        x = rng.standard_normal(50)
        y = rng.standard_normal(50)
        assert autocorrelation_constraint(x, y) >= 0.0

    def test_short_series_returns_zero(self):
        # len < 3 → autocorr function returns 0.0 for both; violation = max(0, |0-0| - tol) = 0
        assert autocorrelation_constraint(np.array([1.0, 2.0]), np.array([1.0, 2.0])) == 0.0


# ---------------------------------------------------------------------------
# annual_stats_constraint
# ---------------------------------------------------------------------------

class TestAnnualStatsConstraint:
    def test_identical_distributions_no_violation(self):
        rng = np.random.default_rng(4)
        x = rng.lognormal(mean=5.0, sigma=0.5, size=120)
        assert annual_stats_constraint(x, x) == pytest.approx(0.0)

    def test_doubled_mean_produces_violation(self):
        base = np.ones(120) * 100.0
        doubled = np.ones(120) * 200.0
        # mean_ratio = 2.0, |2.0 - 1.0| = 1.0 > default mean_tolerance=0.5
        violation = annual_stats_constraint(doubled, base, mean_tolerance=0.5)
        assert violation > 0.0

    def test_within_tolerance_no_violation(self):
        base = np.ones(120) * 100.0
        close = np.ones(120) * 130.0  # ratio 1.3, within tolerance=0.5
        assert annual_stats_constraint(close, base, mean_tolerance=0.5) == pytest.approx(0.0)

    def test_violation_is_non_negative(self):
        rng = np.random.default_rng(5)
        x = rng.lognormal(5, 0.5, 120)
        y = rng.lognormal(6, 1.0, 120)
        assert annual_stats_constraint(x, y) >= 0.0

    def test_zero_historical_mean_skips_mean_check(self):
        # h_mean = 0 → mean check is skipped, only CV violation possible
        historical = np.zeros(120)
        synthetic = np.ones(120)
        # Should not raise; cv_violation only
        violation = annual_stats_constraint(synthetic, historical)
        assert violation >= 0.0


# ---------------------------------------------------------------------------
# seasonal_cycle_constraint
# ---------------------------------------------------------------------------

class TestSeasonalCycleConstraint:
    def _make_2d(self, n_years: int, base_value: float = 100.0) -> np.ndarray:
        """Return (n_years, 12) array with given base value."""
        return np.full((n_years, 12), base_value)

    def test_identical_patterns_no_violation(self):
        data = self._make_2d(30)
        assert seasonal_cycle_constraint(data, data, tolerance=0.5) == pytest.approx(0.0)

    def test_doubled_all_months_produces_violation(self):
        hist = self._make_2d(30, base_value=100.0)
        synth = self._make_2d(30, base_value=250.0)  # ratio=2.5, > 1+0.5=1.5
        assert seasonal_cycle_constraint(synth, hist, tolerance=0.5) > 0.0

    def test_within_tolerance_no_violation(self):
        hist = self._make_2d(30, base_value=100.0)
        synth = self._make_2d(30, base_value=120.0)  # ratio=1.2, within tol=0.5
        assert seasonal_cycle_constraint(synth, hist, tolerance=0.5) == pytest.approx(0.0)

    def test_violation_is_non_negative(self):
        rng = np.random.default_rng(6)
        hist = rng.lognormal(5, 0.3, (30, 12))
        synth = rng.lognormal(6, 0.5, (30, 12))
        assert seasonal_cycle_constraint(synth, hist) >= 0.0
