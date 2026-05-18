"""Safety-net tests for src.discovery.satisficing_labels.

Phase 0 regression net for the src/ -> src/discovery/ move. Covers the
SatisficingDefinition direction logic and apply_labels partitioning.
"""

import numpy as np
import pandas as pd
import pytest

from src.discovery.satisficing_labels import (
    SatisficingDefinition,
    apply_labels,
    VALID_DIRECTIONS,
)


class TestSatisficingDefinitionApply:
    def test_directions_partition_at_threshold(self):
        s = pd.Series([0.9, 1.0, 1.1])
        le = SatisficingDefinition("d", "m", 1.0, "le").apply(s)
        lt = SatisficingDefinition("d", "m", 1.0, "lt").apply(s)
        ge = SatisficingDefinition("d", "m", 1.0, "ge").apply(s)
        gt = SatisficingDefinition("d", "m", 1.0, "gt").apply(s)
        assert le.tolist() == [1, 1, 0]   # boundary included
        assert lt.tolist() == [1, 0, 0]   # boundary excluded
        assert ge.tolist() == [0, 1, 1]
        assert gt.tolist() == [0, 0, 1]
        # le and gt are complementary; lt and ge are complementary.
        assert (le.astype(int) + gt.astype(int)).tolist() == [1, 1, 1]
        assert (lt.astype(int) + ge.astype(int)).tolist() == [1, 1, 1]

    def test_invalid_direction_raises(self):
        with pytest.raises(ValueError):
            SatisficingDefinition("d", "m", 1.0, "between").apply(pd.Series([1.0]))

    def test_valid_directions_constant(self):
        assert VALID_DIRECTIONS == {"le", "lt", "ge", "gt"}


class TestApplyLabels:
    def _bank(self):
        return pd.DataFrame(
            {"realization_id": ["0", "1", "2"], "m1": [0.1, 0.5, np.nan]}
        ).set_index("realization_id")

    def test_long_form_schema_and_one_label_per_row(self):
        defs = [SatisficingDefinition("d1", "m1", 0.3, "le")]
        out = apply_labels(self._bank(), defs)
        assert list(out.columns) == [
            "definition_id", "realization_id", "metric",
            "metric_value", "threshold", "direction", "y",
        ]
        # NaN metric row dropped; 2 valid rows; y is strictly {0, 1}.
        assert len(out) == 2
        assert set(out["y"].unique()).issubset({0, 1})
        assert out["y"].tolist() == [1, 0]

    def test_missing_metric_warns_and_skips(self):
        defs = [SatisficingDefinition("d1", "nope", 0.3, "le")]
        with pytest.warns(UserWarning):
            out = apply_labels(self._bank(), defs)
        assert out.empty

    def test_empty_manifest_returns_empty_with_schema(self):
        out = apply_labels(self._bank(), [])
        assert out.empty
        assert "definition_id" in out.columns
