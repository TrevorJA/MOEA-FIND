"""Library-and-subsample baseline for MOEA-FIND comparison.

Implements the conventional approach to structured drought ensemble generation:
  1. Generate a large library of synthetic traces via SynHydro Kirsch
  2. Characterize each trace by drought metrics (SSI-based)
  3. Subsample using space-filling designs (LHS, Sobol) in drought space

This serves as the primary comparison baseline for MOEA-FIND (DD-09).
The key question: does MOEA-FIND's direct optimization produce better
coverage than post-hoc subsampling from a large library?

Usage:
    from src.hydrology.library import LibraryGenerator

    lib = LibraryGenerator(kirsch_gen, n_years_out=15, ssi_timescale=3)
    lib.generate_library(n_traces=10000, seed=42)
    lib.characterize(metric_set="primary")
    subset_lhs = lib.subsample(method="lhs", n_select=100)
    subset_sobol = lib.subsample(method="sobol", n_select=100)
"""

from typing import Dict, List, Literal, Optional, Tuple

import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class LibraryGenerator:
    """Generate and subsample a large library of synthetic traces.

    Wraps SynHydro's KirschGenerator for bulk trace generation, then
    computes drought characteristics and supports space-filling subsampling.

    Attributes:
        kirsch_gen: Fitted SynHydro KirschGenerator instance.
        n_years_out: Years per synthetic trace.
        ssi_timescale: SSI accumulation period (1, 3, 6, 12).
        traces: Generated traces, shape (n_traces, n_years_out * 12).
        characteristics: DataFrame of drought metrics per trace.
    """

    def __init__(
        self,
        kirsch_gen,
        n_years_out: int = 15,
        ssi_timescale: int = 3,
        reference_flows: Optional[np.ndarray] = None,
    ):
        """Initialize LibraryGenerator.

        Args:
            kirsch_gen: Fitted SynHydro KirschGenerator instance.
            n_years_out: Length of each synthetic trace in years.
            ssi_timescale: SSI accumulation period for drought characterization.
            reference_flows: 1D array of historical monthly flows for SSI fitting.
                If None, must call set_reference_flows() before characterize().
        """
        self.kirsch_gen = kirsch_gen
        self.n_years_out = n_years_out
        self.ssi_timescale = ssi_timescale
        self.reference_flows = reference_flows

        self.traces: Optional[np.ndarray] = None
        self.characteristics: Optional[pd.DataFrame] = None
        self._ssi_calculator = None

    def set_reference_flows(self, reference_flows: np.ndarray) -> None:
        """Set historical reference flows for SSI fitting.

        Args:
            reference_flows: 1D array of historical monthly flows.
        """
        self.reference_flows = reference_flows
        self._ssi_calculator = None

    def generate_library(
        self,
        n_traces: int = 10000,
        seed: Optional[int] = 42,
    ) -> np.ndarray:
        """Generate a library of synthetic traces using SynHydro Kirsch.

        Args:
            n_traces: Number of traces to generate.
            seed: Random seed for reproducibility.

        Returns:
            Array of shape (n_traces, n_years_out * 12) with monthly flows.
        """
        logger.info(
            "Generating %d traces (%d years each)...", n_traces, self.n_years_out
        )

        ensemble = self.kirsch_gen.generate(
            n_realizations=n_traces,
            n_years=self.n_years_out,
            seed=seed,
        )

        # SynHydro returns an Ensemble object with data_by_realization dict
        traces = []
        for real_id in sorted(ensemble.data_by_realization.keys()):
            df = ensemble.data_by_realization[real_id]
            arr = df.values.flatten()
            traces.append(arr[:self.n_years_out * 12])

        self.traces = np.array(traces)
        logger.info("Generated library: shape %s", self.traces.shape)
        return self.traces

    def _get_ssi_calculator(self):
        """Get or create a pre-fitted SSI calculator."""
        if self._ssi_calculator is not None:
            return self._ssi_calculator

        if self.reference_flows is None:
            raise ValueError(
                "reference_flows must be set before characterization. "
                "Call set_reference_flows() or pass to constructor."
            )

        from src.metrics.objectives import make_ssi_calculator, flows_to_series

        calc = make_ssi_calculator(timescale=self.ssi_timescale)
        ref_series = flows_to_series(self.reference_flows, start_date="1950-10-01")
        calc.fit(ref_series)
        self._ssi_calculator = calc
        return calc

    def characterize(
        self,
        metric_set="primary",
    ) -> pd.DataFrame:
        """Compute drought characteristics for all library traces.

        Args:
            metric_set: Either a preset name (e.g. ``"primary"``), a
                single metric name, a sequence of metric names, or a
                tuple of :class:`src.metrics.drought_metrics.DroughtMetric`
                instances. All chars-dict keys are computed regardless;
                this controls which subset is highlighted in the
                summary log and which is treated as the default
                objective set in :meth:`subsample`.

        Returns:
            DataFrame with one row per trace, columns for each metric.
        """
        if self.traces is None:
            raise ValueError("No traces generated. Call generate_library() first.")

        from src.metrics.objectives import (
            flows_to_series,
            compute_ssi_drought_characteristics,
        )
        from src.metrics.drought_metrics import metric_names, resolve_metric_set

        metrics = resolve_metric_set(metric_set)
        names = metric_names(metrics)

        calc = self._get_ssi_calculator()
        n_traces = len(self.traces)

        logger.info("Characterizing %d traces (SSI-%d)...", n_traces, self.ssi_timescale)

        all_chars = []
        for i in range(n_traces):
            syn_series = flows_to_series(
                self.traces[i], start_date="2100-01-01"
            )
            ssi_values = calc.transform(syn_series)
            chars = compute_ssi_drought_characteristics(
                ssi_values, monthly_flows=self.traces[i]
            )
            all_chars.append(chars)

            if (i + 1) % 1000 == 0:
                logger.info("  Characterized %d / %d traces", i + 1, n_traces)

        self.characteristics = pd.DataFrame(all_chars)
        self._metric_set = metrics
        # Stable string view for backwards-compatible callers.
        self._objective_keys = names

        # Log summary for the active metric set.
        for m, name in zip(metrics, names):
            if name in self.characteristics.columns:
                vals = self.characteristics[name]
                logger.info(
                    "  %s: min=%.2f, max=%.2f, mean=%.2f",
                    name, vals.min(), vals.max(), vals.mean(),
                )

        return self.characteristics

    def subsample(
        self,
        method: Literal["lhs", "sobol"] = "lhs",
        n_select: int = 100,
        metric_set=None,
        seed: int = 42,
    ) -> Dict:
        """Subsample the library using a space-filling design in drought space.

        Maps the design points in drought characteristic space to the nearest
        library members. Returns indices, characteristics, and traces.

        Args:
            method: Subsampling method ("lhs" or "sobol").
            n_select: Number of traces to select.
            metric_set: Drought metric set defining the subsampling space.
                Either a preset name, a single metric name, a sequence of
                metric names, or a tuple of :class:`DroughtMetric`. When
                ``None``, falls back to the metric set used in
                :meth:`characterize`.
            seed: Random seed for the space-filling design.

        Returns:
            Dict with keys: indices, characteristics, traces,
            design_points, selected_points, coverage, bounds,
            objective_keys.
        """
        if self.characteristics is None:
            raise ValueError(
                "Library not characterized. Call characterize() first."
            )

        from src.metrics.drought_metrics import metric_names, resolve_metric_set

        if metric_set is None:
            objective_keys = self._objective_keys
        else:
            metrics = resolve_metric_set(metric_set)
            objective_keys = metric_names(metrics)

        from src.discovery.analysis import (
            coverage_metrics,
            generate_lhs_samples,
            generate_sobol_samples,
        )
        from scipy.spatial import KDTree

        # Extract library points in objective space
        lib_points = self.characteristics[list(objective_keys)].values
        lb = lib_points.min(axis=0)
        ub = lib_points.max(axis=0)

        # Ensure non-zero range
        rng = ub - lb
        rng[rng == 0] = 1.0
        ub = lb + rng

        d = len(objective_keys)

        # Generate space-filling design in the library's feasible region
        if method == "lhs":
            design = generate_lhs_samples(n_select, d, lb, ub, seed=seed)
        elif method == "sobol":
            design = generate_sobol_samples(n_select, d, lb, ub, seed=seed)
            design = design[:n_select]
        else:
            raise ValueError(f"Unknown method: {method!r}. Use 'lhs' or 'sobol'.")

        # Find nearest library member for each design point
        tree = KDTree(lib_points)
        _, nearest_idx = tree.query(design)

        # Deduplicate (same library member may be nearest to multiple design points)
        unique_idx = np.unique(nearest_idx)
        logger.info(
            "Subsampled %d traces (%d unique) from %d library via %s",
            n_select, len(unique_idx), len(lib_points), method.upper(),
        )

        selected_chars = self.characteristics.iloc[unique_idx].reset_index(drop=True)
        selected_traces = self.traces[unique_idx]
        selected_points = lib_points[unique_idx]

        # Compute coverage metrics on the selected subset
        cov = coverage_metrics(selected_points, lb, ub)

        return {
            "method": method,
            "n_requested": n_select,
            "n_selected": len(unique_idx),
            "indices": unique_idx,
            "characteristics": selected_chars,
            "traces": selected_traces,
            "design_points": design,
            "selected_points": selected_points,
            "coverage": cov,
            "bounds": {"lb": lb.tolist(), "ub": ub.tolist()},
            "objective_keys": list(objective_keys),
        }

    def save_library(self, path: str) -> None:
        """Save generated library and characteristics to disk.

        Args:
            path: Base path (without extension). Saves:
                {path}_traces.npy - synthetic traces
                {path}_chars.csv - drought characteristics
        """
        if self.traces is None:
            raise ValueError("No traces to save.")

        np.save(f"{path}_traces.npy", self.traces)
        logger.info("Saved traces: %s_traces.npy", path)

        if self.characteristics is not None:
            self.characteristics.to_csv(f"{path}_chars.csv", index=False)
            logger.info("Saved characteristics: %s_chars.csv", path)

    def load_library(self, path: str) -> None:
        """Load a previously saved library from disk.

        Args:
            path: Base path (without extension). Loads:
                {path}_traces.npy - synthetic traces
                {path}_chars.csv - drought characteristics
        """
        self.traces = np.load(f"{path}_traces.npy")
        logger.info("Loaded traces: shape %s", self.traces.shape)

        chars_path = f"{path}_chars.csv"
        try:
            self.characteristics = pd.read_csv(chars_path)
            logger.info("Loaded characteristics: %d rows", len(self.characteristics))
        except FileNotFoundError:
            logger.warning("No characteristics file found at %s", chars_path)
            self.characteristics = None
