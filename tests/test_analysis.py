"""Unit tests for src.analysis — numpy/scipy only, no SynHydro required."""

import numpy as np
import pytest

from src.analysis import coverage_metrics, generate_lhs_samples, generate_sobol_samples


# ---------------------------------------------------------------------------
# coverage_metrics
# ---------------------------------------------------------------------------

class TestCoverageMetrics:
    def test_returns_expected_keys(self):
        points = np.array([[0.0, 0.0], [0.5, 0.5], [1.0, 1.0]])
        lb = np.array([0.0, 0.0])
        ub = np.array([1.0, 1.0])
        m = coverage_metrics(points, lb, ub)
        for key in ("n_points", "dimensions", "L2_star_discrepancy", "nn_mean", "nn_cv"):
            assert key in m, f"missing key: {key}"

    def test_n_points_and_dimensions(self):
        points = np.ones((5, 3))
        m = coverage_metrics(points, np.zeros(3), np.ones(3))
        assert m["n_points"] == 5
        assert m["dimensions"] == 3

    def test_perfect_2d_grid_has_low_nn_cv(self):
        # 4x4 regular grid in [0,1]^2
        xs = np.linspace(0, 1, 4)
        grid = np.array([[x, y] for x in xs for y in xs])
        m = coverage_metrics(grid, np.zeros(2), np.ones(2))
        # Perfect grid should have near-zero NN_CV
        assert m["nn_cv"] < 0.2

    def test_single_point_no_nn_keys(self):
        point = np.array([[0.5, 0.5]])
        m = coverage_metrics(point, np.zeros(2), np.ones(2))
        assert m["n_points"] == 1
        # Nearest-neighbor stats require >1 point
        assert "nn_mean" not in m

    def test_discrepancy_is_positive(self):
        rng = np.random.default_rng(7)
        points = rng.uniform(0, 10, (20, 2))
        m = coverage_metrics(points, np.zeros(2), np.full(2, 10.0))
        assert m["L2_star_discrepancy"] > 0.0

    def test_points_outside_bounds_are_clipped(self):
        # Points that exceed [lb, ub] are clipped to [0,1] before discrepancy
        points = np.array([[-1.0, 0.5], [11.0, 0.5]])
        m = coverage_metrics(points, np.zeros(2), np.ones(2))
        assert "L2_star_discrepancy" in m


# ---------------------------------------------------------------------------
# generate_lhs_samples
# ---------------------------------------------------------------------------

class TestGenerateLHSSamples:
    def test_correct_shape(self):
        samples = generate_lhs_samples(20, 3, np.zeros(3), np.ones(3))
        assert samples.shape == (20, 3)

    def test_values_within_bounds(self):
        lb = np.array([2.0, -5.0])
        ub = np.array([8.0, 0.0])
        samples = generate_lhs_samples(50, 2, lb, ub)
        assert np.all(samples >= lb)
        assert np.all(samples <= ub)

    def test_reproducible_with_seed(self):
        lb, ub = np.zeros(2), np.ones(2)
        s1 = generate_lhs_samples(10, 2, lb, ub, seed=0)
        s2 = generate_lhs_samples(10, 2, lb, ub, seed=0)
        np.testing.assert_array_equal(s1, s2)

    def test_different_seeds_differ(self):
        lb, ub = np.zeros(2), np.ones(2)
        s1 = generate_lhs_samples(10, 2, lb, ub, seed=1)
        s2 = generate_lhs_samples(10, 2, lb, ub, seed=2)
        assert not np.allclose(s1, s2)


# ---------------------------------------------------------------------------
# generate_sobol_samples
# ---------------------------------------------------------------------------

class TestGenerateSobolSamples:
    def test_correct_shape_power_of_2(self):
        # n=8 is a power of 2, shape should be (8, 2)
        samples = generate_sobol_samples(8, 2, np.zeros(2), np.ones(2))
        assert samples.shape == (8, 2)

    def test_values_within_bounds(self):
        lb = np.array([1.0, 10.0])
        ub = np.array([5.0, 20.0])
        samples = generate_sobol_samples(16, 2, lb, ub)
        assert np.all(samples >= lb)
        assert np.all(samples <= ub)

    def test_shape_rounded_to_next_power_of_2(self):
        # n=10 → ceil(log2(10)) = 4 → 2^4 = 16
        samples = generate_sobol_samples(10, 2, np.zeros(2), np.ones(2))
        assert samples.shape == (16, 2)

    def test_reproducible_with_seed(self):
        lb, ub = np.zeros(2), np.ones(2)
        s1 = generate_sobol_samples(8, 2, lb, ub, seed=0)
        s2 = generate_sobol_samples(8, 2, lb, ub, seed=0)
        np.testing.assert_array_equal(s1, s2)
