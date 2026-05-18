"""Safety-net tests for src.hydrology.precompute_trace_series.

Phase 0 regression net for the src/ -> src/hydrology/ move. Builds a
tiny synthetic Pywr-DRB-shaped HDF5 and asserts shape, value range,
idempotence, and the realization-id length guard.
"""

import h5py
import numpy as np
import pandas as pd
import pytest

from src.hydrology.precompute_trace_series import (
    _decode_time,
    compute_annual_nyc_min_storage_frac,
    NYC_STORAGE_CAPACITIES,
    NYC_TOTAL_CAPACITY,
)


def _write_synthetic_hdf5(path, n_days=730, n_traces=4, seed=0):
    rng = np.random.default_rng(seed)
    dates = pd.date_range("2000-01-01", periods=n_days, freq="D")
    with h5py.File(path, "w") as f:
        f.create_dataset(
            "time",
            data=np.array([str(d.date()) for d in dates], dtype="S10"),
        )
        for res, cap in NYC_STORAGE_CAPACITIES.items():
            # Storage strictly within (0, cap] so fraction stays in (0, 1].
            arr = rng.uniform(0.2 * cap, cap, size=(n_days, n_traces))
            f.create_dataset(f"reservoir_{res}", data=arr)
    return dates


class TestDecodeTime:
    def test_bytes_decoded(self):
        with h5py.File("/dev/null", "w", driver="core", backing_store=False) as f:
            ds = f.create_dataset(
                "t", data=np.array(["2001-01-01", "2001-01-02"], dtype="S10")
            )
            idx = _decode_time(ds)
        assert isinstance(idx, pd.DatetimeIndex)
        assert idx[0] == pd.Timestamp("2001-01-01")


class TestComputeAnnualNYCMinStorageFrac:
    def test_shape_range_and_idempotent(self, tmp_path):
        h5 = tmp_path / "pywrdrb_output.hdf5"
        _write_synthetic_hdf5(h5, n_days=730, n_traces=4)
        df = compute_annual_nyc_min_storage_frac(h5)
        # 730 days starting 2000-01-01 spans calendar years 2000 and 2001.
        assert list(df.index) == [2000, 2001]
        assert df.shape == (2, 4)
        assert df.index.name == "year"
        vals = df.to_numpy()
        assert np.all(vals > 0.0) and np.all(vals <= 1.0)
        assert not np.isnan(vals).any()
        # Idempotent: recompute yields identical frame.
        df2 = compute_annual_nyc_min_storage_frac(h5)
        pd.testing.assert_frame_equal(df, df2)

    def test_custom_realization_ids(self, tmp_path):
        h5 = tmp_path / "out.hdf5"
        _write_synthetic_hdf5(h5, n_days=400, n_traces=3)
        df = compute_annual_nyc_min_storage_frac(
            h5, realization_ids=["a", "b", "c"]
        )
        assert list(df.columns) == ["a", "b", "c"]

    def test_realization_id_length_mismatch_raises(self, tmp_path):
        h5 = tmp_path / "out.hdf5"
        _write_synthetic_hdf5(h5, n_days=400, n_traces=3)
        with pytest.raises(ValueError):
            compute_annual_nyc_min_storage_frac(h5, realization_ids=["a", "b"])

    def test_total_capacity_consistent(self):
        assert NYC_TOTAL_CAPACITY == float(sum(NYC_STORAGE_CAPACITIES.values()))
