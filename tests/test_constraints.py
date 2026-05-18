"""Tests for src/constraints.py (ConstraintConfig + two-tier violations).

These tests exercise the Block-B rewrite of the plausibility-constraint module.
They do not require the calibration JSON to exist: fixtures build a synthetic
``ConstraintConfig`` in-memory with tight, well-known tolerances so every
branch of the violation logic is reachable without IO.

The old hand-tolerance functions (autocorrelation_constraint,
annual_stats_constraint, seasonal_cycle_constraint) no longer exist — the
calibrated two-tier violations replace them.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from src.optimization.constraints import (
    ConstraintConfig,
    ConstraintResult,
    Violation,
    annual_cv_violation,
    annual_mean_violation,
    compute_all_constraints,
    lag1_ac_violation,
    seasonal_cycle_violation,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _historical_cycle() -> np.ndarray:
    """Canonical cannonsville-like seasonal cycle (high spring, low late summer)."""
    return np.array([
        300.0, 320.0, 450.0, 700.0, 900.0, 700.0,
        450.0, 350.0, 300.0, 260.0, 250.0, 280.0,
    ])


def _make_trace(
    cycle: np.ndarray,
    n_years: int = 20,
    noise_sd: float = 20.0,
    seed: int = 0,
) -> np.ndarray:
    """Generate a plausible monthly flow trace from a seasonal cycle + noise."""
    rng = np.random.default_rng(seed)
    years = np.tile(cycle, (n_years, 1))
    noise = rng.normal(0.0, noise_sd, size=years.shape)
    return np.clip(years + noise, 1.0, None)


@pytest.fixture
def base_cfg():
    cycle = _historical_cycle()
    hist_2d = _make_trace(cycle, n_years=70, noise_sd=30.0, seed=42)
    hist_1d = hist_2d.flatten()
    hist_annual_totals = hist_2d.sum(axis=1)
    hist_monthly_means = hist_2d.mean(axis=0)

    return ConstraintConfig(
        T_years=20,
        hist_annual_mean=float(hist_annual_totals.mean()),
        hist_annual_cv=float(hist_annual_totals.std(ddof=1) / hist_annual_totals.mean()),
        hist_lag1_ac_monthly=float(
            np.corrcoef(hist_1d[:-1], hist_1d[1:])[0, 1]
        ),
        hist_non_drought_mean=float(hist_1d.mean()),
        hist_monthly_means=hist_monthly_means.tolist(),
        annual_mean_tol=0.15,
        annual_cv_tol=0.25,
        lag1_ac_tol=0.10,
        non_drought_mean_tol=0.15,
        seasonal_cycle_tol=0.20,
        hard_multiplier=2.0,
        soft_weight=1.0,
        enabled={
            "annual_mean": True,
            "annual_cv": True,
            "lag1_ac_monthly": True,
            "non_drought_mean": False,
            "seasonal_cycle": True,
        },
    )


# ---------------------------------------------------------------------------
# ConstraintConfig round-trip
# ---------------------------------------------------------------------------

class TestConstraintConfigIO:

    def test_to_dict_round_trip(self, base_cfg):
        d = base_cfg.to_dict()
        cfg2 = ConstraintConfig.from_dict(d)
        assert cfg2.T_years == base_cfg.T_years
        assert cfg2.annual_mean_tol == pytest.approx(base_cfg.annual_mean_tol)
        assert cfg2.enabled == base_cfg.enabled

    def test_json_round_trip(self, tmp_path: Path, base_cfg):
        path = tmp_path / "cfg.json"
        base_cfg.to_json(path)
        loaded = json.loads(path.read_text())
        cfg2 = ConstraintConfig.from_dict(loaded)
        assert cfg2.hist_monthly_means == base_cfg.hist_monthly_means

    def test_from_calibration_json(self, tmp_path: Path):
        payload = {
            "cannonsville_T20": {
                "T_years": 20,
                "historical": {
                    "annual_mean": 5000.0,
                    "annual_cv": 0.3,
                    "lag1_ac_monthly": 0.85,
                    "non_drought_mean": 450.0,
                    "monthly_means": [400.0] * 12,
                },
                "tolerances": {
                    "annual_mean": 0.14,
                    "annual_cv": 0.28,
                    "lag1_ac_monthly": 0.09,
                    "non_drought_mean": 0.13,
                    "seasonal_cycle_max_dev": 0.19,
                },
            }
        }
        path = tmp_path / "calib.json"
        path.write_text(json.dumps(payload))
        cfg = ConstraintConfig.from_calibration_json(path, T_years=20)
        assert cfg.annual_mean_tol == pytest.approx(0.14)
        assert cfg.lag1_ac_tol == pytest.approx(0.09)
        assert cfg.seasonal_cycle_tol == pytest.approx(0.19)
        assert len(cfg.hist_monthly_means) == 12


# ---------------------------------------------------------------------------
# Two-tier violation logic — individual constraints
# ---------------------------------------------------------------------------

class TestAnnualMeanViolation:

    def test_clean_trace_is_feasible(self, base_cfg):
        cycle = _historical_cycle()
        trace = _make_trace(cycle, n_years=20, noise_sd=20.0, seed=1).flatten()
        v = annual_mean_violation(trace, base_cfg)
        assert v.hard_violation == 0.0
        assert v.soft_penalty == 0.0
        assert v.name == "annual_mean"

    def test_halved_trace_is_hard_infeasible(self, base_cfg):
        cycle = _historical_cycle()
        trace = 0.4 * _make_trace(cycle, n_years=20, noise_sd=10.0, seed=2).flatten()
        v = annual_mean_violation(trace, base_cfg)
        assert v.hard_violation > 0.0
        assert v.deviation > base_cfg.annual_mean_tol * base_cfg.hard_multiplier

    def test_moderate_deviation_is_soft_only(self, base_cfg):
        cycle = _historical_cycle()
        trace = 1.20 * _make_trace(cycle, n_years=20, noise_sd=5.0, seed=3).flatten()
        v = annual_mean_violation(trace, base_cfg)
        assert v.hard_violation == 0.0
        assert 0.0 < v.soft_penalty <= 1.0

    def test_soft_penalty_monotonic(self, base_cfg):
        cycle = _historical_cycle()
        base_trace = _make_trace(cycle, n_years=20, noise_sd=1.0, seed=4).flatten()
        scales = [1.16, 1.22, 1.28]
        penalties = [
            annual_mean_violation(scale * base_trace, base_cfg).soft_penalty
            for scale in scales
        ]
        assert penalties[0] < penalties[1] < penalties[2]


class TestAnnualCvViolation:

    def test_clean_trace_is_feasible(self, base_cfg):
        cycle = _historical_cycle()
        trace = _make_trace(cycle, n_years=20, noise_sd=30.0, seed=5).flatten()
        v = annual_cv_violation(trace, base_cfg)
        assert v.hard_violation == 0.0

    def test_inflated_variability_flagged(self, base_cfg):
        cycle = _historical_cycle()
        years = _make_trace(cycle, n_years=20, noise_sd=5.0, seed=6)
        amp = np.array([0.3 if i % 2 == 0 else 1.7 for i in range(20)])
        trace = (years * amp[:, None]).flatten()
        v = annual_cv_violation(trace, base_cfg)
        assert v.hard_violation > 0.0 or v.soft_penalty > 0.0


class TestLag1AcViolation:

    def test_clean_trace_is_feasible(self, base_cfg):
        cycle = _historical_cycle()
        trace = _make_trace(cycle, n_years=20, noise_sd=25.0, seed=7).flatten()
        v = lag1_ac_violation(trace, base_cfg)
        assert v.hard_violation == 0.0

    def test_iid_noise_flagged(self, base_cfg):
        rng = np.random.default_rng(8)
        trace = rng.uniform(100.0, 900.0, size=240)
        v = lag1_ac_violation(trace, base_cfg)
        assert v.hard_violation > 0.0


class TestSeasonalCycleViolation:

    def test_clean_trace_is_feasible(self, base_cfg):
        cycle = _historical_cycle()
        trace = _make_trace(cycle, n_years=20, noise_sd=15.0, seed=9)
        v = seasonal_cycle_violation(trace, base_cfg)
        assert v.hard_violation == 0.0

    def test_flat_trace_is_flagged(self, base_cfg):
        trace = np.full((20, 12), 450.0)
        v = seasonal_cycle_violation(trace, base_cfg)
        assert v.hard_violation > 0.0


# ---------------------------------------------------------------------------
# Aggregated constraint result
# ---------------------------------------------------------------------------

class TestComputeAllConstraints:

    def test_clean_trace_all_feasible(self, base_cfg):
        cycle = _historical_cycle()
        trace_2d = _make_trace(cycle, n_years=20, noise_sd=20.0, seed=10)
        trace_1d = trace_2d.flatten()
        result = compute_all_constraints(trace_1d, trace_2d, base_cfg)
        assert isinstance(result, ConstraintResult)
        assert result.feasible is True
        assert all(h == 0.0 for h in result.hard_violations)
        assert result.soft_penalty_total == pytest.approx(0.0, abs=1e-9)

    def test_degenerate_trace_flagged_on_multiple(self, base_cfg):
        rng = np.random.default_rng(11)
        one_year = rng.uniform(50.0, 100.0, size=12)
        trace_2d = np.tile(one_year, (20, 1))
        trace_1d = trace_2d.flatten()
        result = compute_all_constraints(trace_1d, trace_2d, base_cfg)
        assert result.feasible is False
        assert len(result.hard_violations) == 4
        assert sum(result.hard_violations) > 0.0

    def test_non_drought_requires_ssi_calc(self, base_cfg):
        cfg = ConstraintConfig.from_dict(base_cfg.to_dict())
        cfg.enabled["non_drought_mean"] = True
        cycle = _historical_cycle()
        trace_2d = _make_trace(cycle, n_years=20, noise_sd=20.0, seed=12)
        with pytest.raises(ValueError, match="non_drought_mean"):
            compute_all_constraints(trace_2d.flatten(), trace_2d, cfg)

    def test_soft_weight_applied(self, base_cfg):
        base_cfg.soft_weight = 2.5
        cycle = _historical_cycle()
        trace_2d = 1.20 * _make_trace(cycle, n_years=20, noise_sd=2.0, seed=13)
        trace_1d = trace_2d.flatten()
        result = compute_all_constraints(trace_1d, trace_2d, base_cfg)
        assert result.soft_penalty_weighted == pytest.approx(
            2.5 * result.soft_penalty_total
        )
