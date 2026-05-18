"""Unit tests for src.metrics.drought_metrics — registry, presets, anti-ideal placement."""

import numpy as np
import pytest

from src.metrics.drought_metrics import (
    AntiIdealRule,
    DroughtMetric,
    PRESETS,
    REGISTRY,
    compute_anti_ideal,
    metric_epsilons,
    metric_labels,
    metric_names,
    resolve_metric_set,
)


# ---------------------------------------------------------------------------
# Registry shape
# ---------------------------------------------------------------------------


class TestRegistry:
    def test_primary_metrics_are_continuous_and_present(self):
        for name in PRESETS["primary"]:
            assert name in REGISTRY, f"{name} missing from REGISTRY"
            assert not REGISTRY[name].is_cyclic, (
                f"primary preset must contain only non-cyclic metrics; "
                f"{name} is cyclic"
            )

    def test_primary_excludes_clustered_metrics(self):
        # User direction: the primary preset must not contain duration
        # or peak month (both cluster at integer/discrete values).
        primary = set(PRESETS["primary"])
        assert "mean_duration" not in primary
        assert "peak_severity_month" not in primary

    def test_legacy_preset_holds_old_defaults(self):
        # Reproducibility of pre-2026-04-27 runs requires the legacy
        # preset to remain available with the original three keys.
        assert PRESETS["legacy"] == (
            "mean_duration",
            "mean_avg_severity",
            "peak_severity_month",
        )

    def test_every_preset_resolves(self):
        for name in PRESETS:
            metrics = resolve_metric_set(name)
            assert all(isinstance(m, DroughtMetric) for m in metrics)

    def test_metric_names_are_unique_within_each_preset(self):
        for name, names in PRESETS.items():
            assert len(names) == len(set(names)), f"duplicate in preset {name}"

    def test_each_metric_has_well_formed_metadata(self):
        for name, m in REGISTRY.items():
            assert m.name == name
            assert m.label, f"{name} missing label"
            assert m.units, f"{name} missing units"
            assert m.epsilon > 0, f"{name} has non-positive epsilon"
            assert isinstance(m.anti_ideal_rule, AntiIdealRule)
            assert callable(m.extract)


# ---------------------------------------------------------------------------
# resolve_metric_set
# ---------------------------------------------------------------------------


class TestResolveMetricSet:
    def test_preset_name_returns_tuple_of_metrics(self):
        metrics = resolve_metric_set("primary")
        assert isinstance(metrics, tuple)
        assert len(metrics) == len(PRESETS["primary"])
        assert metric_names(metrics) == PRESETS["primary"]

    def test_single_metric_name(self):
        metrics = resolve_metric_set("mean_severity")
        assert len(metrics) == 1
        assert metrics[0].name == "mean_severity"

    def test_sequence_of_names(self):
        metrics = resolve_metric_set(["mean_severity", "mean_magnitude"])
        assert metric_names(metrics) == ("mean_severity", "mean_magnitude")

    def test_unknown_name_raises(self):
        with pytest.raises(KeyError):
            resolve_metric_set("not_a_metric")

    def test_unknown_in_sequence_raises(self):
        with pytest.raises(KeyError):
            resolve_metric_set(["mean_severity", "not_a_metric"])

    def test_pass_through_metric_instances(self):
        # A tuple of DroughtMetric instances is returned as-is.
        existing = resolve_metric_set("primary")
        again = resolve_metric_set(existing)
        assert again == existing

    def test_empty_sequence_raises(self):
        with pytest.raises(KeyError):
            resolve_metric_set([])


# ---------------------------------------------------------------------------
# compute_anti_ideal
# ---------------------------------------------------------------------------


class TestComputeAntiIdeal:
    def test_constant_rule_uses_metric_constant(self):
        # time_in_drought_fraction has CONSTANT rule with default 1.0.
        m = REGISTRY["time_in_drought_fraction"]
        assert m.anti_ideal_rule is AntiIdealRule.CONSTANT
        out = compute_anti_ideal((m,), hist_chars={}, headroom=1.5)
        assert out.shape == (1,)
        assert out[0] == pytest.approx(1.0)

    def test_cyclic_rule_uses_period_times_headroom(self):
        m = REGISTRY["peak_severity_month"]
        assert m.anti_ideal_rule is AntiIdealRule.CYCLIC_HEADROOM
        out = compute_anti_ideal((m,), hist_chars={}, headroom=1.5)
        assert out[0] == pytest.approx(18.0)

    def test_headroom_times_max_uses_max_partner(self):
        # mean_severity has max_partner = "worst_severity".
        m = REGISTRY["mean_severity"]
        out = compute_anti_ideal(
            (m,),
            hist_chars={"mean_severity": 1.0, "worst_severity": 2.5},
            headroom=2.0,
        )
        assert out[0] == pytest.approx(5.0)

    def test_headroom_times_max_falls_back_to_self_when_no_partner(self):
        # frequency has no max_partner declared; it uses its own value.
        m = REGISTRY["frequency"]
        out = compute_anti_ideal(
            (m,),
            hist_chars={"frequency": 4.0},
            headroom=1.5,
        )
        assert out[0] == pytest.approx(6.0)

    def test_feasible_maxes_override_for_headroom_metrics(self):
        m = REGISTRY["mean_severity"]
        out = compute_anti_ideal(
            (m,),
            hist_chars={"worst_severity": 2.0},
            headroom=1.5,
            feasible_maxes={"mean_severity": 3.0},
        )
        assert out[0] == pytest.approx(4.5)

    def test_feasible_maxes_ignored_for_cyclic_and_constant(self):
        cyclic = REGISTRY["peak_severity_month"]
        constant = REGISTRY["time_in_drought_fraction"]
        out = compute_anti_ideal(
            (cyclic, constant),
            hist_chars={},
            headroom=1.5,
            feasible_maxes={"peak_severity_month": 99.0,
                            "time_in_drought_fraction": 99.0},
        )
        assert out[0] == pytest.approx(18.0)
        assert out[1] == pytest.approx(1.0)

    def test_zero_max_falls_back_to_default(self):
        # All-zero historical case (no events) should yield a non-zero D*.
        m = REGISTRY["mean_severity"]
        out = compute_anti_ideal(
            (m,),
            hist_chars={"worst_severity": 0.0},
            headroom=1.5,
        )
        assert out[0] > 0.0

    def test_assembles_per_metric_into_vector(self):
        metric_set = resolve_metric_set("primary")
        out = compute_anti_ideal(
            metric_set,
            hist_chars={
                "worst_severity": 2.0,
                "max_magnitude": 8.0,
            },
            headroom=1.5,
        )
        assert out.shape == (3,)
        # Each entry is non-negative; the time_in_drought entry equals 1.0.
        assert all(out >= 0)
        assert out[2] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_metric_names_round_trips(self):
        metrics = resolve_metric_set("extreme_event")
        assert metric_names(metrics) == PRESETS["extreme_event"]

    def test_metric_labels_returns_strings(self):
        labels = metric_labels(resolve_metric_set("primary"))
        assert all(isinstance(l, str) and l for l in labels)

    def test_metric_epsilons_match_registry(self):
        metrics = resolve_metric_set("primary")
        eps = metric_epsilons(metrics)
        assert eps == tuple(m.epsilon for m in metrics)


# ---------------------------------------------------------------------------
# Sign convention: more drought ⇒ larger objective (DD-11)
# ---------------------------------------------------------------------------


class TestSignConvention:
    def test_q10_flow_extractor_returns_negated_quantile(self):
        # The q10_flow metric is keyed off "q10_flow_neg" so that larger
        # objective corresponds to lower flow (more drought stress).
        m = REGISTRY["q10_flow"]
        chars = {"q10_flow_neg": -150.0}
        # extractor returns the value as-is from the chars dict.
        assert m.extract(chars) == pytest.approx(-150.0)

    def test_extractor_returns_zero_for_missing_key(self):
        m = REGISTRY["mean_severity"]
        assert m.extract({}) == 0.0
