"""Unit tests for KirschBorgWrapper — no SynHydro required.

A minimal mock KirschGenerator is constructed from scratch so that
KirschBorgWrapper's __init__, property, and mapping logic can be tested
without any SynHydro import.
"""

import numpy as np
import pytest

from src.hydrology.kirsch_wrapper import KirschBorgWrapper


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_gen(
    n_historic_years: int = 50,
    n_sites: int = 1,
) -> object:
    """Return a minimal mock KirschGenerator-like object.

    Provides the attributes that KirschBorgWrapper.__init__ reads:
      Z_h, mean_month, std_month, n_historic_years, n_sites
    """

    class MockKirschGenerator:
        pass

    gen = MockKirschGenerator()
    gen.n_historic_years = n_historic_years
    gen.n_sites = n_sites
    # Z_h: (n_years, 12, n_sites) standardized residuals
    rng = np.random.default_rng(42)
    gen.Z_h = rng.standard_normal((n_historic_years, 12, n_sites))
    gen.mean_month = rng.uniform(3.0, 7.0, (12, n_sites))
    gen.std_month = rng.uniform(0.1, 0.5, (12, n_sites))
    return gen


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------

class TestKirschBorgWrapperInit:
    def test_valid_index_mode(self):
        gen = _make_mock_gen()
        wrapper = KirschBorgWrapper(gen, mode="index", n_years_out=30)
        assert wrapper.mode == "index"
        assert wrapper.n_years_out == 30

    def test_valid_residual_mode(self):
        gen = _make_mock_gen()
        wrapper = KirschBorgWrapper(gen, mode="residual", n_years_out=20)
        assert wrapper.mode == "residual"
        assert wrapper.sorted_residuals is not None

    def test_invalid_mode_raises_value_error(self):
        gen = _make_mock_gen()
        with pytest.raises(ValueError, match="mode must be"):
            KirschBorgWrapper(gen, mode="bad_mode")

    def test_unfitted_generator_missing_Z_h_raises(self):
        class BareGen:
            n_historic_years = 50
            n_sites = 1
            mean_month = np.ones((12, 1))
            std_month = np.ones((12, 1))

        with pytest.raises(AttributeError, match="fitted before wrapping"):
            KirschBorgWrapper(BareGen(), mode="index")

    def test_unfitted_generator_missing_mean_month_raises(self):
        class BareGen:
            n_historic_years = 50
            n_sites = 1
            Z_h = np.ones((50, 12, 1))
            std_month = np.ones((12, 1))

        with pytest.raises(AttributeError, match="mean_month"):
            KirschBorgWrapper(BareGen(), mode="index")

    def test_attributes_stored_correctly(self):
        gen = _make_mock_gen(n_historic_years=40, n_sites=2)
        wrapper = KirschBorgWrapper(gen, mode="index", n_years_out=25)
        assert wrapper.n_years_hist == 40
        assert wrapper.n_sites == 2
        assert wrapper.kirsch_gen is gen


# ---------------------------------------------------------------------------
# n_dvs property
# ---------------------------------------------------------------------------

class TestNDvs:
    def test_index_mode_n_dvs(self):
        gen = _make_mock_gen()
        wrapper = KirschBorgWrapper(gen, mode="index", n_years_out=30)
        # index mode: (n_years_out + 1) * 12
        assert wrapper.n_dvs == 31 * 12

    def test_residual_mode_n_dvs(self):
        gen = _make_mock_gen()
        wrapper = KirschBorgWrapper(gen, mode="residual", n_years_out=30)
        # residual mode: n_years_out * 12
        assert wrapper.n_dvs == 30 * 12

    def test_n_dvs_scales_with_n_years_out(self):
        gen = _make_mock_gen()
        for n in (10, 20, 50):
            w_idx = KirschBorgWrapper(gen, mode="index", n_years_out=n)
            w_res = KirschBorgWrapper(gen, mode="residual", n_years_out=n)
            assert w_idx.n_dvs == (n + 1) * 12
            assert w_res.n_dvs == n * 12


# ---------------------------------------------------------------------------
# _map_indices
# ---------------------------------------------------------------------------

class TestMapIndices:
    def setup_method(self):
        self.gen = _make_mock_gen(n_historic_years=50)
        self.wrapper = KirschBorgWrapper(self.gen, mode="index", n_years_out=30)

    def test_zero_maps_to_zero_index(self):
        dvs = np.zeros(self.wrapper.n_dvs)
        indices = self.wrapper._map_indices(dvs)
        assert np.all(indices == 0)

    def test_one_maps_to_last_valid_index(self):
        dvs = np.ones(self.wrapper.n_dvs)
        indices = self.wrapper._map_indices(dvs)
        # floor(1.0 * 50) = 50, clipped to 49
        assert np.all(indices == 49)

    def test_all_indices_within_valid_range(self):
        rng = np.random.default_rng(10)
        dvs = rng.uniform(0, 1, self.wrapper.n_dvs)
        indices = self.wrapper._map_indices(dvs)
        assert np.all(indices >= 0)
        assert np.all(indices < self.gen.n_historic_years)

    def test_output_dtype_is_int(self):
        dvs = np.full(self.wrapper.n_dvs, 0.5)
        indices = self.wrapper._map_indices(dvs)
        assert indices.dtype in (np.int32, np.int64, int)

    def test_midpoint_maps_to_middle_index(self):
        dvs = np.full(self.wrapper.n_dvs, 0.5)
        indices = self.wrapper._map_indices(dvs)
        # floor(0.5 * 50) = 25
        assert np.all(indices == 25)


# ---------------------------------------------------------------------------
# sorted_residuals shape (residual mode)
# ---------------------------------------------------------------------------

class TestSortedResiduals:
    def test_sorted_residuals_shape(self):
        gen = _make_mock_gen(n_historic_years=50, n_sites=2)
        wrapper = KirschBorgWrapper(gen, mode="residual", n_years_out=20)
        # sorted along axis=0 of Z_h, shape unchanged
        assert wrapper.sorted_residuals.shape == (50, 12, 2)

    def test_sorted_residuals_are_sorted_along_axis0(self):
        gen = _make_mock_gen(n_historic_years=30, n_sites=1)
        wrapper = KirschBorgWrapper(gen, mode="residual", n_years_out=10)
        sr = wrapper.sorted_residuals
        for m in range(12):
            col = sr[:, m, 0]
            assert np.all(col[:-1] <= col[1:]), f"Not sorted for month {m}"
