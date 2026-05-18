"""Tests for the Phase-6 metric-family dispatcher (src.metrics.dispatch)."""

from src.metrics.dispatch import (
    resolve_metric_family,
    chars_fn_for,
    SHORT_BLOCK_MAX_YEARS,
)


class TestResolveMetricFamily:
    def test_short_block_for_t1_t2(self):
        for t in (1, 2):
            fam, handle = resolve_metric_family(t)
            assert fam == "short_block"
            assert hasattr(handle, "compute_short_block_metrics")

    def test_candidates_mode_for_long_t(self):
        fam, handle = resolve_metric_family(3, mode="candidates")
        assert fam == "extended"
        assert hasattr(handle, "compute_all_candidates")

    def test_registry_default_for_long_t(self):
        fam, handle = resolve_metric_family(3)
        assert fam == "registry"
        # REGISTRY is the production drought-metric mapping.
        assert len(handle) > 0

    def test_threshold_constant(self):
        assert SHORT_BLOCK_MAX_YEARS == 2
        assert resolve_metric_family(SHORT_BLOCK_MAX_YEARS)[0] == "short_block"
        assert resolve_metric_family(SHORT_BLOCK_MAX_YEARS + 1)[0] == "registry"


class TestCharsFnFor:
    def test_long_t_returns_none(self):
        # Longer traces use the inline pre-fitted-SSI path in the runner.
        assert chars_fn_for(10, None) is None

    def test_short_t_returns_callable(self):
        import numpy as np
        rng = np.random.default_rng(0)
        monthly = rng.lognormal(6.0, 0.3, size=40 * 12)
        fn = chars_fn_for(1, monthly)
        assert callable(fn)
