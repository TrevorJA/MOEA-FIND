"""Safety-net tests for src.hydrology.historical_blocks (pure resampling, no SynHydro).

Added in Phase 0 of the spring-cleaning reorg as a regression net so the
subsequent module move (src/ -> src/hydrology/) cannot silently break
historical-block resampling. Asserts shape, determinism, and the
subset-membership invariant.
"""

import numpy as np
import pytest

from src.hydrology.historical_blocks import (
    resample_historical_blocks,
    resample_disjoint_tilings,
    resample_historical_blocks_2d,
)


def _record(n_years: int) -> np.ndarray:
    # Deterministic, distinct per-month values so membership is checkable.
    return np.arange(n_years * 12, dtype=float)


class TestResampleHistoricalBlocks:
    def test_block_shape_and_count_stride1(self):
        rec = _record(10)
        blocks = resample_historical_blocks(rec, T_years=3, stride=1)
        # floor((10 - 3) / 1) + 1 = 8 overlapping blocks
        assert len(blocks) == 8
        assert all(b.shape == (36,) for b in blocks)

    def test_stride_reduces_count(self):
        rec = _record(10)
        n1 = len(resample_historical_blocks(rec, T_years=3, stride=1))
        n3 = len(resample_historical_blocks(rec, T_years=3, stride=3))
        assert n3 < n1

    def test_subset_membership_and_chronological(self):
        rec = _record(8)
        blocks = resample_historical_blocks(rec, T_years=2, stride=1)
        recset = set(rec.tolist())
        for b in blocks:
            assert set(b.tolist()).issubset(recset)
        # First block is the earliest window.
        np.testing.assert_array_equal(blocks[0], rec[:24])

    def test_determinism(self):
        rec = _record(12)
        a = resample_historical_blocks(rec, T_years=4, stride=2)
        b = resample_historical_blocks(rec, T_years=4, stride=2)
        assert len(a) == len(b)
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)

    def test_record_shorter_than_T_raises(self):
        with pytest.raises(ValueError):
            resample_historical_blocks(_record(2), T_years=5)

    def test_bad_stride_raises(self):
        with pytest.raises(ValueError):
            resample_historical_blocks(_record(5), T_years=2, stride=0)

    def test_non_1d_raises(self):
        with pytest.raises(ValueError):
            resample_historical_blocks(_record(5).reshape(5, 12), T_years=2)


class TestResampleDisjointTilings:
    def test_within_tiling_no_year_repeats(self):
        rec = _record(10)
        tilings = resample_disjoint_tilings(rec, T_years=3)
        assert len(tilings) == 3  # default n_offsets == T_years
        for tiling in tilings:
            seen = []
            for blk in tiling:
                seen.extend(blk.tolist())
            assert len(seen) == len(set(seen))  # disjoint within a tiling

    def test_short_record_raises(self):
        with pytest.raises(ValueError):
            resample_disjoint_tilings(_record(2), T_years=5)


class TestResampleHistoricalBlocks2d:
    def test_2d_block_shape(self):
        rec2d = _record(10).reshape(10, 12)
        blocks = resample_historical_blocks_2d(rec2d, T_years=3, stride=1)
        assert len(blocks) == 8
        assert all(b.shape == (3, 12) for b in blocks)

    def test_bad_shape_raises(self):
        with pytest.raises(ValueError):
            resample_historical_blocks_2d(_record(10).reshape(20, 6), T_years=2)
