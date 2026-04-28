"""Unit tests for src.sensitivity — Ishigami-function validation.

The Ishigami function (Ishigami & Homma, 1990) is the standard SA test:

    f(X) = sin(X1) + a sin(X2)^2 + b X3^4 sin(X1),  X ~ U[-pi, pi]^3
    a = 7,  b = 0.1

Analytic first-order Sobol indices are known: roughly
    S1_1 ≈ 0.314,  S1_2 ≈ 0.442,  S1_3 ≈ 0.0
so the importance ranking is X2 > X1 > X3. Every SA method must reproduce
that ranking on a sufficiently large sample.

The tests use a 5000-row LHS sample so SALib's index estimators have
room to converge. ``num_resamples`` / ``n_bootstrap`` are kept small to
keep the suite fast.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("SALib")

from src.sensitivity import (  # noqa: E402
    METHODS,
    apply_method_selection_criterion,
    bootstrap_indices,
    bootstrap_rank_stability,
    compute_delta,
    compute_pawn,
    compute_rbd_fast,
    convergence_curve,
    cross_method_rank_corr,
    cross_outcome_rank_corr,
    drop_low_variance_factors,
)


# ---------------------------------------------------------------------------
# Sample builder
# ---------------------------------------------------------------------------


def _ishigami_sample(n: int = 5000, seed: int = 0):
    """Latin-hypercube sample on [-pi, pi]^3 evaluated through Ishigami."""
    rng = np.random.default_rng(seed)
    # Standard LHS: each column is a jittered permutation of [0, n-1]/n
    perms = np.column_stack([rng.permutation(n) for _ in range(3)])
    u = (perms + rng.random((n, 3))) / n  # values in [0, 1)
    X = -np.pi + 2.0 * np.pi * u
    a, b = 7.0, 0.1
    Y = np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2 + b * X[:, 2] ** 4 * np.sin(X[:, 0])
    return X, Y, ["x1", "x2", "x3"]


def _top_k(df: pd.DataFrame, headline_col: str) -> str:
    """Return the factor name with the largest headline-index value."""
    return str(df[headline_col].idxmax())


# ---------------------------------------------------------------------------
# Per-method ranking on Ishigami
# ---------------------------------------------------------------------------


class TestIshigamiRanking:
    """Each SA method must rank Ishigami factors X2 > X1 > X3."""

    def test_delta(self):
        X, Y, names = _ishigami_sample(n=4000, seed=0)
        df = compute_delta(X, Y, names, num_resamples=20, seed=0)
        assert _top_k(df, "delta") == "x2"
        # X3 has zero first-order contribution; it should be ranked last
        # by the moment-independent index too.
        assert df["delta"].idxmin() == "x3"

    def test_pawn(self):
        X, Y, names = _ishigami_sample(n=4000, seed=0)
        df = compute_pawn(X, Y, names, S=10, seed=0, n_bootstrap=0)
        assert _top_k(df, "median") == "x2"
        assert df["median"].idxmin() == "x3"

    def test_rbd_fast(self):
        X, Y, names = _ishigami_sample(n=4000, seed=0)
        df = compute_rbd_fast(X, Y, names, M=10, num_resamples=20, seed=0)
        assert _top_k(df, "S1") == "x2"
        # Analytic S1_3 = 0; allow numerical noise.
        assert abs(float(df.loc["x3", "S1"])) < 0.1


# ---------------------------------------------------------------------------
# Cross-method rank correlation must agree on Ishigami
# ---------------------------------------------------------------------------


def test_cross_method_rank_corr_high_on_ishigami():
    X, Y, names = _ishigami_sample(n=4000, seed=0)
    results = {
        "delta": compute_delta(X, Y, names, num_resamples=20, seed=0),
        "pawn": compute_pawn(X, Y, names, S=10, seed=0, n_bootstrap=0),
        "rbd_fast": compute_rbd_fast(X, Y, names, M=10, num_resamples=20, seed=0),
    }
    rho = cross_method_rank_corr(results)
    # Every off-diagonal Spearman should be ≥ 0.5 on a clean Ishigami
    # 4000-sample (in practice all near 1.0 because the ranking is
    # identical).
    off_diag = rho.where(~np.eye(len(rho), dtype=bool))
    assert (off_diag.dropna(how="all").stack() >= 0.5).all()
    # Self-correlations are exactly 1.0
    np.testing.assert_allclose(np.diag(rho.values), 1.0)


# ---------------------------------------------------------------------------
# Cross-outcome rank correlation works as advertised
# ---------------------------------------------------------------------------


def test_cross_outcome_rank_corr_identical_outcomes():
    """Two copies of the same outcome must yield rho = 1."""
    X, Y, names = _ishigami_sample(n=2000, seed=0)
    df = compute_delta(X, Y, names, num_resamples=10, seed=0)
    per_outcome = {"a": df, "b": df.copy()}
    rho = cross_outcome_rank_corr(per_outcome, method="delta")
    np.testing.assert_allclose(rho.values, 1.0)


# ---------------------------------------------------------------------------
# Bootstrap indices: shape and CI ordering
# ---------------------------------------------------------------------------


def test_bootstrap_indices_returns_ordered_ci():
    X, Y, names = _ishigami_sample(n=1500, seed=1)
    method = METHODS["delta"]
    df = bootstrap_indices(method, X, Y, names,
                           n_bootstrap=20, seed=0,
                           method_kwargs={"num_resamples": 10})
    assert list(df.columns) == ["mean", "ci_lo", "ci_hi"]
    assert (df["ci_lo"] <= df["ci_hi"] + 1e-12).all()


def test_bootstrap_rank_stability_attrs_present():
    X, Y, names = _ishigami_sample(n=1500, seed=2)
    method = METHODS["delta"]
    df = bootstrap_rank_stability(method, X, Y, names,
                                  n_bootstrap=20, seed=0,
                                  method_kwargs={"num_resamples": 10})
    for k in ("rank_spearman_median",
              "rank_spearman_iqr_lo",
              "rank_spearman_iqr_hi"):
        assert k in df.attrs
    # On Ishigami the ranking is stable; median Spearman should be high.
    assert df.attrs["rank_spearman_median"] > 0.5


# ---------------------------------------------------------------------------
# Convergence curve shape
# ---------------------------------------------------------------------------


def test_convergence_curve_shape():
    X, Y, names = _ishigami_sample(n=2000, seed=3)
    method = METHODS["delta"]
    sizes = [200, 500]
    df = convergence_curve(method, X, Y, names,
                           sizes=sizes, n_replicates=4, seed=0,
                           method_kwargs={"num_resamples": 10})
    expected_n = (len(sizes) + 1) * len(names)  # +1 for full-sample bucket
    assert len(df) == expected_n
    assert set(df.columns) >= {"factor", "n", "mean", "p05", "p95", "replicates"}
    # Every factor appears at every size
    for f in names:
        assert (df["factor"] == f).sum() == len(sizes) + 1


# ---------------------------------------------------------------------------
# Low-variance factor dropping
# ---------------------------------------------------------------------------


def test_drop_low_variance_factors_drops_constant_column():
    X = pd.DataFrame({
        "varies": np.random.default_rng(0).normal(0.0, 1.0, size=200),
        "constant": np.full(200, 3.14),
    })
    kept, dropped = drop_low_variance_factors(X, cv_threshold=0.05)
    assert "varies" in kept.columns
    assert "constant" in dropped


def test_drop_low_variance_factors_keeps_all_when_threshold_zero():
    X = pd.DataFrame({
        "x1": np.random.default_rng(0).normal(0.0, 1.0, size=100),
        "x2": np.random.default_rng(1).normal(2.0, 0.05, size=100),
    })
    kept, dropped = drop_low_variance_factors(X, cv_threshold=0.0)
    assert dropped == []
    assert list(kept.columns) == ["x1", "x2"]


# ---------------------------------------------------------------------------
# Selection criterion: Ishigami passes; constructed-failure case fails
# ---------------------------------------------------------------------------


def test_selection_criterion_passes_on_ishigami():
    X, Y, names = _ishigami_sample(n=3000, seed=4)
    method = METHODS["delta"]
    boot = bootstrap_indices(method, X, Y, names,
                             n_bootstrap=30, seed=0,
                             method_kwargs={"num_resamples": 10})
    rank = bootstrap_rank_stability(method, X, Y, names,
                                    n_bootstrap=30, seed=0,
                                    method_kwargs={"num_resamples": 10})
    # Synthetic cross-method comparison: pretend PAWN agrees perfectly.
    cross_rho = pd.Series({"pawn": 1.0, "rbd_fast": 1.0})
    decision = apply_method_selection_criterion(
        bootstrap_df=boot, rank_stability_df=rank,
        cross_method_rho=cross_rho,
    )
    assert decision["passes"] is True
    assert decision["top_factor"] == "x2"
    assert decision["failing_conditions"] == []


def test_selection_criterion_fails_when_cross_method_disagrees():
    X, Y, names = _ishigami_sample(n=2000, seed=5)
    method = METHODS["delta"]
    boot = bootstrap_indices(method, X, Y, names,
                             n_bootstrap=20, seed=0,
                             method_kwargs={"num_resamples": 10})
    rank = bootstrap_rank_stability(method, X, Y, names,
                                    n_bootstrap=20, seed=0,
                                    method_kwargs={"num_resamples": 10})
    # Cross-method correlation forced below threshold.
    cross_rho = pd.Series({"pawn": 0.1, "rbd_fast": 0.0})
    decision = apply_method_selection_criterion(
        bootstrap_df=boot, rank_stability_df=rank,
        cross_method_rho=cross_rho,
        cross_method_threshold=0.7,
    )
    assert decision["passes"] is False
    assert "cross_method_disagreement" in decision["failing_conditions"]
