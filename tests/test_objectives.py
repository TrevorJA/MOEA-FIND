"""Unit tests for src.objectives — pure math, no SynHydro required."""

import numpy as np
import pytest

from src.drought_metrics import REGISTRY, resolve_metric_set
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
    """drought_objectives implements the DD-11 L1 Device:
        f_j(x)     = D_j(x)                     for j = 1..K (raw metric)
        f_{K+1}(x) = ||D(x) - D*||_1            (Manhattan distance)

    The first K entries are raw drought characteristics (minimised toward 0
    by Borg). The (K+1)th entry is the Manhattan distance to the anti-ideal
    D*. Under D_j <= D*_j (enforced by compute_ssi_anti_ideal placing D*
    outside the feasible region), sum_{j=1}^{K+1} f_j is constant across
    feasible solutions, so every feasible point is Pareto non-dominated and
    Borg's epsilon-box archive tiles the drought-characteristic space.

    This is the load-bearing formulation of MOEA-FIND — see
    manuscript/design_decisions.md §DD-11. Do NOT re-introduce a
    ``|D_j - D*_j|`` variant of the first K objectives; that collapses the
    Pareto front to the single point nearest D*.
    """

    def _chars(self, duration=3.0, avg_severity=1.5):
        return {"mean_duration": duration, "mean_avg_severity": avg_severity}

    def test_output_shape_is_k_plus_1(self):
        objs = drought_objectives(
            self._chars(),
            anti_ideal=np.array([5.0, 2.0]),
        )
        assert objs.shape == (3,)

    def test_first_k_entries_are_raw_metrics(self):
        # DD-11: f_j = D_j for j = 1..K (NOT |D_j - D*_j|).
        anti_ideal = np.array([10.0, 3.0])
        objs = drought_objectives(self._chars(4.0, 1.0), anti_ideal)
        assert objs[0] == pytest.approx(4.0)
        assert objs[1] == pytest.approx(1.0)

    def test_last_entry_is_manhattan_distance_to_anti_ideal(self):
        anti_ideal = np.array([5.0, 2.0])
        objs = drought_objectives(self._chars(3.0, 1.5), anti_ideal)
        expected = manhattan_norm(np.array([3.0, 1.5]), anti_ideal)
        assert objs[-1] == pytest.approx(expected)
        # Explicit check: |3 - 5| + |1.5 - 2| = 2.5.
        assert objs[-1] == pytest.approx(2.5)

    def test_hyperplane_identity_holds_inside_feasible_region(self):
        # DD-11: sum_{j=1..K+1} f_j == sum_j D*_j  whenever D_j <= D*_j.
        anti_ideal = np.array([10.0, 3.0])
        for duration, severity in [(0.0, 0.0), (4.0, 1.0), (9.0, 2.5), (10.0, 3.0)]:
            objs = drought_objectives(
                self._chars(duration, severity), anti_ideal
            )
            assert np.sum(objs) == pytest.approx(np.sum(anti_ideal)), (
                f"hyperplane identity broken at D=({duration}, {severity}): "
                f"sum(objs)={np.sum(objs)}, sum(D*)={np.sum(anti_ideal)}"
            )

    def test_at_anti_ideal_last_objective_is_zero(self):
        # When D_j = D*_j for all j, the Manhattan distance is zero and all
        # of the "drought magnitude" lives in the first K raw objectives.
        anti_ideal = np.array([8.0, 4.0])
        objs = drought_objectives(self._chars(8.0, 4.0), anti_ideal)
        assert objs[0] == pytest.approx(8.0)
        assert objs[1] == pytest.approx(4.0)
        assert objs[-1] == pytest.approx(0.0)

    def test_no_drought_corner_norm_equals_sum_of_anti_ideal(self):
        # At the opposite corner, D_j = 0 for all j, so f_{K+1} = sum(D*).
        anti_ideal = np.array([5.0, 2.0])
        objs = drought_objectives(self._chars(0.0, 0.0), anti_ideal)
        assert objs[:-1] == pytest.approx(np.zeros(2))
        assert objs[-1] == pytest.approx(np.sum(anti_ideal))

    def test_missing_key_defaults_to_zero_metric(self):
        # Missing key → D_j = 0, so f_j = 0 and f_{K+1} picks up the full
        # anti-ideal contribution for that dimension.
        anti_ideal = np.array([5.0, 2.0])
        objs = drought_objectives({}, anti_ideal=anti_ideal)
        assert objs[0] == pytest.approx(0.0)
        assert objs[1] == pytest.approx(0.0)
        assert objs[-1] == pytest.approx(np.sum(anti_ideal))

    def test_anti_ideal_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="anti_ideal must be length"):
            drought_objectives(
                self._chars(), anti_ideal=np.array([1.0]),
                objective_keys=("mean_duration", "mean_avg_severity"),
            )

    def test_cyclic_month_anti_ideal_outside_calendar(self):
        # peak_severity_month ∈ [0, 12]; D*_month = 18 (set by
        # compute_anti_ideal as 12 × 1.5) keeps the hyperplane
        # assumption D_j <= D*_j intact even for cyclic axes.
        chars = {"mean_duration": 5.0, "peak_severity_month": 11.0}
        anti_ideal = np.array([10.0, 18.0])
        objs = drought_objectives(
            chars, anti_ideal=anti_ideal,
            objective_keys=("mean_duration", "peak_severity_month"),
        )
        assert objs[0] == pytest.approx(5.0)
        assert objs[1] == pytest.approx(11.0)
        # Manhattan: |5 - 10| + |11 - 18| = 5 + 7 = 12.
        assert objs[-1] == pytest.approx(12.0)
        # Hyperplane identity still holds with cyclic-month D*.
        assert np.sum(objs) == pytest.approx(np.sum(anti_ideal))

    def test_metric_set_passes_through_drought_metric_instances(self):
        # New API: drought_objectives can take a tuple of DroughtMetric
        # instances directly. Each metric's ``extract`` pulls the value
        # out of the chars dict, so the same algebraic identities hold.
        metrics = resolve_metric_set(["mean_severity", "mean_magnitude"])
        chars = {"mean_severity": 1.5, "mean_magnitude": 4.0}
        anti_ideal = np.array([3.0, 8.0])
        objs = drought_objectives(chars, anti_ideal, objective_keys=metrics)
        assert objs[0] == pytest.approx(1.5)
        assert objs[1] == pytest.approx(4.0)
        assert objs[-1] == pytest.approx(5.5)  # |1.5-3| + |4-8|
        # DD-11 hyperplane identity holds for the metric-set form too.
        assert np.sum(objs) == pytest.approx(np.sum(anti_ideal))


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
