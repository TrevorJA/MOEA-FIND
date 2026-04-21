"""Tests for src/constraints_dv.py (DV-uniformity constraint, ablation arm).

These tests pin down the behaviour the ablation experiment relies on:

    - both statistics (L2-star discrepancy, KS) are bounded and behave
      monotonically on known-extreme inputs;
    - ``compute_dv_constraint`` emits a ``ConstraintResult`` of the shape
      ``run_experiment`` expects when ``n_constrs = 1``;
    - the two-tier soft/hard band machinery inherited from
      ``src.constraints._two_tier`` is wired correctly — soft penalty 0 inside
      tolerance, hard_violation > 0 past the hard cut-off.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.constraints import ConstraintResult, _two_tier
from src.constraints_dv import (
    DVUniformityConfig,
    compute_dv_constraint,
    dv_anderson_darling_statistic,
    dv_ks_statistic,
    dv_l2_star_statistic,
    statistic_fn,
)


# ---------------------------------------------------------------------------
# Statistic functions
# ---------------------------------------------------------------------------

class TestStatistics:
    def test_uniform_rng_draw_scores_small_on_both(self):
        rng = np.random.default_rng(20260420)
        dvs = rng.uniform(0.0, 1.0, size=240)
        l2 = dv_l2_star_statistic(dvs)
        ks = dv_ks_statistic(dvs)
        # N=240 U[0,1] draw: L2-star typically < 1e-2, KS < 0.1.
        assert 0.0 <= l2 < 5e-2, f"L2-star too large for uniform draw: {l2}"
        assert 0.0 <= ks < 0.2, f"KS too large for uniform draw: {ks}"

    def test_constant_half_vector_is_degenerate(self):
        dvs = np.full(240, 0.5)
        l2 = dv_l2_star_statistic(dvs)
        ks = dv_ks_statistic(dvs)
        # Constant at 0.5 => empirical CDF is a step function at 0.5;
        # KS distance from uniform CDF is exactly 0.5.
        assert ks == pytest.approx(0.5, abs=1e-6)
        # L2-star is strictly larger for a degenerate point set than for a
        # well-spread one — use a looser lower bound than the uniform draw.
        assert l2 > 1e-2

    def test_clustered_near_zero_inflates_both(self):
        # Drought-like DVs: concentrate in the low-quantile tail.
        rng = np.random.default_rng(0)
        dvs = rng.uniform(0.0, 0.2, size=240)
        l2 = dv_l2_star_statistic(dvs)
        ks = dv_ks_statistic(dvs)
        # Compared to the uniform case above, both should be much larger.
        rng2 = np.random.default_rng(1)
        dvs_uni = rng2.uniform(0.0, 1.0, size=240)
        l2_uni = dv_l2_star_statistic(dvs_uni)
        ks_uni = dv_ks_statistic(dvs_uni)
        assert l2 > l2_uni
        assert ks > ks_uni
        # KS bounded in [0, 1].
        assert 0.0 <= ks <= 1.0

    def test_statistics_bounded(self):
        rng = np.random.default_rng(7)
        for _ in range(20):
            dvs = rng.uniform(0.0, 1.0, size=rng.integers(30, 300))
            assert 0.0 <= dv_l2_star_statistic(dvs) <= 1.0
            assert 0.0 <= dv_ks_statistic(dvs) <= 1.0

    def test_statistic_fn_dispatch(self):
        dvs = np.linspace(0.0, 1.0, 100)
        assert statistic_fn("l2_star")(dvs) == dv_l2_star_statistic(dvs)
        assert statistic_fn("ks")(dvs) == dv_ks_statistic(dvs)
        assert (
            statistic_fn("ad")(dvs) == dv_anderson_darling_statistic(dvs)
        )
        with pytest.raises(ValueError):
            statistic_fn("hurst")

    def test_empty_input_is_zero(self):
        empty = np.array([], dtype=float)
        assert dv_l2_star_statistic(empty) == 0.0
        assert dv_ks_statistic(empty) == 0.0
        assert dv_anderson_darling_statistic(empty) == 0.0

    def test_ad_uniform_draw_is_small(self):
        rng = np.random.default_rng(20260420)
        dvs = rng.uniform(0.0, 1.0, size=240)
        a2 = dv_anderson_darling_statistic(dvs)
        # Under the null (true U[0,1] draw) A² ≈ 1 on average; 95th
        # percentile ≈ 2.492 for any n. Single draw is unlikely to exceed 5.
        assert 0.0 <= a2 < 5.0, f"AD too large for uniform draw: {a2}"

    def test_ad_increases_when_tail_mass_added(self):
        # Monotonicity: starting from a uniform draw and adding tail
        # mass should strictly increase A². This is the behavioural
        # property the constraint relies on.
        rng = np.random.default_rng(7)
        baseline = rng.uniform(0.0, 1.0, size=240)
        tail_heavy = baseline.copy()
        tail_heavy[:30] = rng.uniform(0.0, 0.02, size=30)
        assert (
            dv_anderson_darling_statistic(tail_heavy)
            > dv_anderson_darling_statistic(baseline)
        )

    def test_ad_constant_vector_is_huge(self):
        # Pathological degenerate input: all DVs at 0.5. A² should be
        # much larger than any uniform-draw value.
        dvs = np.full(240, 0.5)
        a2 = dv_anderson_darling_statistic(dvs)
        assert a2 > 50.0, f"AD on constant 0.5 should be large, got {a2}"


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

class TestConfig:
    def test_rejects_invalid_statistic(self):
        with pytest.raises(ValueError):
            DVUniformityConfig(T_years=20, n_dvs=240,
                               statistic="bogus", tolerance=0.01)

    def test_rejects_negative_tolerance(self):
        with pytest.raises(ValueError):
            DVUniformityConfig(T_years=20, n_dvs=240,
                               statistic="l2_star", tolerance=-1.0)

    def test_rejects_nonpositive_n_dvs(self):
        with pytest.raises(ValueError):
            DVUniformityConfig(T_years=20, n_dvs=0,
                               statistic="l2_star", tolerance=0.01)

    def test_roundtrip_dict(self):
        cfg = DVUniformityConfig(
            T_years=20, n_dvs=240, statistic="ks", tolerance=0.12,
            hard_multiplier=3.0, soft_weight=0.25, enabled=False,
        )
        d = cfg.to_dict()
        back = DVUniformityConfig.from_dict(d)
        assert back.T_years == 20
        assert back.n_dvs == 240
        assert back.statistic == "ks"
        assert back.tolerance == pytest.approx(0.12)
        assert back.hard_multiplier == 3.0
        assert back.soft_weight == 0.25
        assert back.enabled is False


# ---------------------------------------------------------------------------
# compute_dv_constraint
# ---------------------------------------------------------------------------

class TestComputeDVConstraint:
    @staticmethod
    def _cfg(tolerance: float = 0.05, statistic: str = "l2_star",
             hard_multiplier: float = 2.0, soft_weight: float = 0.5,
             enabled: bool = True, n_dvs: int = 240) -> DVUniformityConfig:
        return DVUniformityConfig(
            T_years=20, n_dvs=n_dvs, statistic=statistic,
            tolerance=tolerance, hard_multiplier=hard_multiplier,
            soft_weight=soft_weight, enabled=enabled,
        )

    def test_result_shape(self):
        rng = np.random.default_rng(0)
        dvs = rng.uniform(0.0, 1.0, size=240)
        cfg = self._cfg()
        res = compute_dv_constraint(dvs, cfg)
        assert isinstance(res, ConstraintResult)
        assert len(res.hard_violations) == 1
        assert len(res.violations) == 1
        assert isinstance(res.feasible, bool)
        # run_experiment expects soft_penalty_weighted to be float-compatible.
        assert isinstance(res.soft_penalty_weighted, float)

    def test_disabled_is_feasible(self):
        dvs = np.full(240, 0.5)  # would normally fail hard
        cfg = self._cfg(enabled=False)
        res = compute_dv_constraint(dvs, cfg)
        assert res.feasible
        assert res.hard_violations == [0.0]
        assert res.soft_penalty_total == 0.0
        assert res.soft_penalty_weighted == 0.0

    def test_uniform_draw_is_feasible_with_reasonable_tolerance(self):
        rng = np.random.default_rng(20260420)
        dvs = rng.uniform(0.0, 1.0, size=240)
        cfg = self._cfg(tolerance=0.05)  # generous for N=240
        res = compute_dv_constraint(dvs, cfg)
        assert res.feasible
        assert res.hard_violations == [0.0]
        assert res.soft_penalty_weighted == 0.0

    def test_constant_vector_hard_violates(self):
        # Degenerate DVs => statistic is well above a tight tolerance.
        dvs = np.full(240, 0.5)
        cfg = self._cfg(tolerance=1e-4, hard_multiplier=2.0, statistic="ks")
        res = compute_dv_constraint(dvs, cfg)
        assert not res.feasible
        assert res.hard_violations[0] > 0.0
        # KS = 0.5, hard cut = 2e-4, so hard_violation ≈ 0.4998.
        assert res.hard_violations[0] == pytest.approx(0.5 - 2e-4, abs=1e-3)

    def test_soft_penalty_in_band(self):
        # Pick DVs whose statistic lands between tol and hard_cut.
        # We construct a cfg with large tolerance around a known KS value.
        dvs = np.full(240, 0.5)  # KS exactly 0.5
        cfg = self._cfg(tolerance=0.4, hard_multiplier=2.0, statistic="ks")
        # stat=0.5, tol=0.4, hard_cut=0.8 -> inside soft band.
        res = compute_dv_constraint(dvs, cfg)
        assert res.feasible  # hard_violation == 0
        assert res.hard_violations[0] == 0.0
        assert 0.0 < res.soft_penalty_total <= 1.0
        assert res.soft_penalty_weighted == pytest.approx(
            0.5 * res.soft_penalty_total
        )

    def test_dv_length_mismatch_raises(self):
        cfg = self._cfg(n_dvs=240)
        with pytest.raises(ValueError, match="DV length mismatch"):
            compute_dv_constraint(np.zeros(100), cfg)

    def test_two_tier_quadratic_progression(self):
        # Sanity: _two_tier semantics match the docstring.
        # dev at tol => soft_penalty = 0; at hard_cut => soft_penalty = 1.
        v_at_tol = _two_tier("x", 0.10, tolerance=0.10, hard_multiplier=2.0)
        v_at_hard = _two_tier("x", 0.20, tolerance=0.10, hard_multiplier=2.0)
        v_mid = _two_tier("x", 0.15, tolerance=0.10, hard_multiplier=2.0)
        assert v_at_tol.soft_penalty == 0.0
        assert v_at_hard.soft_penalty == pytest.approx(1.0)
        # Quadratic halfway through the band: (0.5)^2 = 0.25.
        assert v_mid.soft_penalty == pytest.approx(0.25)
