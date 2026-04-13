"""Unit tests for src.objectives — pure math, no SynHydro required."""

import numpy as np
import pytest

from src.objectives import (
    analytic_objectives,
    compute_drought_characteristics,
    compute_drought_events,
    drought_objectives,
    manhattan_norm,
    normalize_to_unit_cube,
)


# ---------------------------------------------------------------------------
# manhattan_norm
# ---------------------------------------------------------------------------

class TestManhattanNorm:
    def test_known_values(self):
        assert manhattan_norm(np.array([1, 2]), np.array([5, 5])) == pytest.approx(7.0)

    def test_zero_distance(self):
        assert manhattan_norm(np.array([3, 3]), np.array([3, 3])) == pytest.approx(0.0)

    def test_1d(self):
        assert manhattan_norm(np.array([0]), np.array([5])) == pytest.approx(5.0)

    def test_negative_coordinates(self):
        assert manhattan_norm(np.array([-1, -1]), np.array([1, 1])) == pytest.approx(4.0)

    def test_returns_float(self):
        result = manhattan_norm(np.array([1.0, 2.0]), np.array([3.0, 4.0]))
        assert isinstance(result, float)


# ---------------------------------------------------------------------------
# analytic_objectives
# ---------------------------------------------------------------------------

class TestAnalyticObjectives:
    def test_2d_case(self):
        objs = analytic_objectives(np.array([1.0, 2.0]), np.array([5.0, 5.0]))
        assert objs.tolist() == pytest.approx([1.0, 2.0, 7.0])

    def test_last_element_equals_manhattan_norm(self):
        dvs = np.array([2.0, 3.0])
        anti_ideal = np.array([6.0, 7.0])
        objs = analytic_objectives(dvs, anti_ideal)
        expected_norm = manhattan_norm(dvs, anti_ideal)
        assert objs[-1] == pytest.approx(expected_norm)

    def test_3d_case(self):
        objs = analytic_objectives(np.array([0.0, 0.0, 0.0]), np.array([3.0, 3.0, 3.0]))
        assert objs.tolist() == pytest.approx([0.0, 0.0, 0.0, 9.0])

    def test_output_length_is_k_plus_1(self):
        for k in (1, 2, 3, 5):
            dvs = np.zeros(k)
            anti_ideal = np.ones(k)
            assert len(analytic_objectives(dvs, anti_ideal)) == k + 1

    def test_hyperplane_property(self):
        """sum of all k+1 objectives == sum(anti_ideal) when dvs <= anti_ideal."""
        anti_ideal = np.array([4.0, 6.0, 8.0])
        # Any point with all components <= anti_ideal lies on the hyperplane
        dvs = np.array([1.0, 3.0, 5.0])
        objs = analytic_objectives(dvs, anti_ideal)
        assert np.sum(objs) == pytest.approx(np.sum(anti_ideal))


# ---------------------------------------------------------------------------
# drought_objectives
# ---------------------------------------------------------------------------

class TestDroughtObjectives:
    def _chars(self, duration=3.0, avg_severity=1.5):
        return {"mean_duration": duration, "mean_avg_severity": avg_severity}

    def test_output_shape_is_k_plus_1(self):
        objs = drought_objectives(
            self._chars(),
            anti_ideal=np.array([5.0, 2.0]),
        )
        assert objs.shape == (3,)

    def test_manhattan_norm_matches(self):
        anti_ideal = np.array([5.0, 2.0])
        objs = drought_objectives(self._chars(3.0, 1.5), anti_ideal)
        expected = manhattan_norm(np.array([3.0, 1.5]), anti_ideal)
        assert objs[-1] == pytest.approx(expected)

    def test_zero_events_norm_equals_sum_of_anti_ideal(self):
        anti_ideal = np.array([5.0, 2.0])
        objs = drought_objectives(self._chars(0.0, 0.0), anti_ideal)
        assert objs[-1] == pytest.approx(np.sum(anti_ideal))

    def test_metric_values_in_first_k_elements(self):
        anti_ideal = np.array([10.0, 3.0])
        objs = drought_objectives(self._chars(4.0, 1.0), anti_ideal)
        assert objs[0] == pytest.approx(4.0)
        assert objs[1] == pytest.approx(1.0)

    def test_missing_key_defaults_to_zero(self):
        chars = {}  # no keys at all
        objs = drought_objectives(chars, anti_ideal=np.array([5.0, 2.0]))
        assert objs[0] == pytest.approx(0.0)
        assert objs[1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# normalize_to_unit_cube
# ---------------------------------------------------------------------------

class TestNormalizeToUnitCube:
    def test_bounds_map_to_0_and_1(self):
        lb = np.array([0.0, 0.0])
        ub = np.array([10.0, 10.0])
        points = np.array([[0.0, 0.0], [10.0, 10.0]])
        normed = normalize_to_unit_cube(points, lb, ub)
        assert normed[0].tolist() == pytest.approx([0.0, 0.0])
        assert normed[1].tolist() == pytest.approx([1.0, 1.0])

    def test_midpoint_maps_to_half(self):
        lb = np.array([0.0, 0.0])
        ub = np.array([10.0, 10.0])
        points = np.array([[5.0, 5.0]])
        normed = normalize_to_unit_cube(points, lb, ub)
        assert normed[0].tolist() == pytest.approx([0.5, 0.5])

    def test_asymmetric_bounds(self):
        lb = np.array([2.0, -4.0])
        ub = np.array([6.0, 0.0])
        points = np.array([[4.0, -2.0]])
        normed = normalize_to_unit_cube(points, lb, ub)
        assert normed[0].tolist() == pytest.approx([0.5, 0.5])


# ---------------------------------------------------------------------------
# compute_drought_events
# ---------------------------------------------------------------------------

class TestComputeDroughtEvents:
    def test_two_events(self):
        flows = np.array([100.0, 50.0, 40.0, 60.0, 100.0, 30.0, 100.0])
        events = compute_drought_events(flows, threshold=55.0)
        assert len(events) == 2

    def test_first_event_duration_and_deficit(self):
        flows = np.array([100.0, 50.0, 40.0, 60.0, 100.0, 30.0, 100.0])
        events = compute_drought_events(flows, threshold=55.0)
        e = events[0]
        assert e["start"] == 1
        assert e["duration"] == 2
        assert e["deficit"] == pytest.approx(20.0)  # (55-50) + (55-40)
        assert e["intensity"] == pytest.approx(10.0)

    def test_no_drought(self):
        flows = np.array([100.0, 200.0, 300.0])
        assert compute_drought_events(flows, threshold=50.0) == []

    def test_all_drought_is_one_event(self):
        flows = np.array([10.0, 20.0, 30.0])
        events = compute_drought_events(flows, threshold=50.0)
        assert len(events) == 1
        assert events[0]["duration"] == 3

    def test_event_at_end_of_series(self):
        flows = np.array([100.0, 20.0, 30.0])
        events = compute_drought_events(flows, threshold=50.0)
        assert len(events) == 1
        e = events[0]
        assert e["start"] == 1
        assert e["duration"] == 2
        assert e["deficit"] == pytest.approx(50.0)  # (50-20)+(50-30)

    def test_single_month_event(self):
        flows = np.array([100.0, 10.0, 100.0])
        events = compute_drought_events(flows, threshold=50.0)
        assert len(events) == 1
        assert events[0]["duration"] == 1


# ---------------------------------------------------------------------------
# compute_drought_characteristics
# ---------------------------------------------------------------------------

class TestComputeDroughtCharacteristics:
    def test_frequency_calculation(self):
        # 2 events over 10 years = 2 events/decade
        flows = np.array([100.0, 50.0, 40.0, 60.0, 100.0, 30.0, 100.0])
        chars = compute_drought_characteristics(flows, threshold=55.0, n_years=10.0)
        assert chars["frequency"] == pytest.approx(2.0)

    def test_mean_duration(self):
        flows = np.array([100.0, 50.0, 40.0, 60.0, 100.0, 30.0, 100.0])
        chars = compute_drought_characteristics(flows, threshold=55.0, n_years=10.0)
        # event 1: duration=2, event 2: duration=1 -> mean=1.5
        assert chars["mean_duration"] == pytest.approx(1.5)

    def test_mean_intensity(self):
        flows = np.array([100.0, 50.0, 40.0, 60.0, 100.0, 30.0, 100.0])
        chars = compute_drought_characteristics(flows, threshold=55.0, n_years=10.0)
        # event 1: intensity=10.0, event 2: intensity=25.0 -> mean=17.5
        assert chars["mean_intensity"] == pytest.approx(17.5)

    def test_no_events_returns_zeros(self):
        flows = np.array([100.0, 200.0, 300.0])
        chars = compute_drought_characteristics(flows, threshold=50.0, n_years=10.0)
        assert chars["frequency"] == 0.0
        assert chars["mean_duration"] == 0.0
        assert chars["n_events"] == 0
