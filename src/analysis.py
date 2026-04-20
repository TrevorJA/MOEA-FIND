"""Post-processing and visualization for MOEA-FIND experiments."""

from typing import Dict, Optional, Sequence, Tuple
import numpy as np


def pick_space_filling_subset(
    points: np.ndarray,
    n: int,
    seed: int = 42,
) -> np.ndarray:
    """Pick ``n`` indices from *points* that are approximately space-filling.

    Used to select a small, reproducible subset of the Pareto archive for
    the dev-ensemble smoke test. Algorithm: draw a Latin Hypercube over
    the axis-aligned bounding box of *points*, then for each LHS anchor
    pick the nearest not-yet-used point from *points*.

    Args:
        points: ``(N, D)`` array of candidate coordinates.
        n: Number of indices to return. Must satisfy ``n <= N``.
        seed: LHS random seed.

    Returns:
        1D array of ``n`` unique indices into *points*.
    """
    from scipy.spatial import cKDTree
    from scipy.stats.qmc import LatinHypercube

    pts = np.asarray(points, dtype=float)
    N, D = pts.shape
    if n > N:
        raise ValueError(f"requested {n} points but only {N} available")
    if n == N:
        return np.arange(N)

    lb = pts.min(axis=0)
    ub = pts.max(axis=0)
    # Guard against degenerate axes (all-equal)
    ub = np.where(ub <= lb, lb + 1e-9, ub)

    sampler = LatinHypercube(d=D, seed=seed)
    anchors = lb + sampler.random(n=n) * (ub - lb)

    tree = cKDTree(pts)
    chosen: list = []
    used = set()
    # For each anchor, query a small k, pick the first not-yet-used hit.
    k_pool = min(max(8, n // 4), N)
    for anchor in anchors:
        _, idxs = tree.query(anchor, k=k_pool)
        idxs = np.atleast_1d(idxs)
        for i in idxs:
            i = int(i)
            if i not in used:
                used.add(i)
                chosen.append(i)
                break
        else:
            # Pool exhausted of unused neighbors; fall back to global search
            mask = np.ones(N, dtype=bool)
            mask[list(used)] = False
            remaining_idx = np.where(mask)[0]
            if len(remaining_idx) == 0:
                break
            d2 = ((pts[remaining_idx] - anchor) ** 2).sum(axis=1)
            pick = int(remaining_idx[np.argmin(d2)])
            used.add(pick)
            chosen.append(pick)
    return np.array(chosen, dtype=int)


def generate_lhs_samples(
    n: int,
    d: int,
    lb: np.ndarray,
    ub: np.ndarray,
    seed: int = 42,
) -> np.ndarray:
    """Generate Latin Hypercube samples in [lb, ub]^d.

    Args:
        n: Number of samples.
        d: Dimensionality.
        lb: Lower bounds (d,).
        ub: Upper bounds (d,).
        seed: Random seed.

    Returns:
        Array of shape (n, d).
    """
    from scipy.stats.qmc import LatinHypercube
    sampler = LatinHypercube(d=d, seed=seed)
    unit_samples = sampler.random(n=n)
    return lb + unit_samples * (ub - lb)


def generate_sobol_samples(
    n: int,
    d: int,
    lb: np.ndarray,
    ub: np.ndarray,
    seed: int = 42,
) -> np.ndarray:
    """Generate Sobol quasi-random samples in [lb, ub]^d.

    Args:
        n: Number of samples (will be rounded to next power of 2).
        d: Dimensionality.
        lb: Lower bounds (d,).
        ub: Upper bounds (d,).
        seed: Random seed.

    Returns:
        Array of shape (n, d).
    """
    from scipy.stats.qmc import Sobol
    sampler = Sobol(d=d, scramble=True, seed=seed)
    # Sobol requires power of 2
    m = int(np.ceil(np.log2(n)))
    unit_samples = sampler.random_base2(m=m)
    return lb + unit_samples * (ub - lb)


def coverage_metrics(
    points: np.ndarray,
    lb: np.ndarray,
    ub: np.ndarray,
) -> Dict[str, float]:
    """Compute coverage quality metrics for a point set.

    Args:
        points: Array of shape (n, d) in original coordinates.
        lb: Lower bounds for normalization.
        ub: Upper bounds for normalization.

    Returns:
        Dict with L2-star discrepancy, centered discrepancy, and
        nearest-neighbor statistics.
    """
    from scipy.stats.qmc import discrepancy
    from scipy.spatial import KDTree

    # Normalize to [0,1]^d
    normed = (points - lb) / (ub - lb)
    normed = np.clip(normed, 0, 1)

    metrics = {}
    metrics["n_points"] = len(points)
    metrics["dimensions"] = points.shape[1]
    metrics["L2_star_discrepancy"] = float(discrepancy(normed, method="L2-star"))

    # Nearest-neighbor distance statistics
    if len(points) > 1:
        tree = KDTree(normed)
        dists, _ = tree.query(normed, k=2)  # k=2 because first neighbor is self
        nn_dists = dists[:, 1]
        metrics["nn_mean"] = float(np.mean(nn_dists))
        metrics["nn_std"] = float(np.std(nn_dists))
        metrics["nn_min"] = float(np.min(nn_dists))
        metrics["nn_max"] = float(np.max(nn_dists))
        # Coefficient of variation: lower = more uniform spacing
        metrics["nn_cv"] = float(np.std(nn_dists) / np.mean(nn_dists))

    return metrics
