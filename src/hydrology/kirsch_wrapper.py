"""Borg-compatible wrapper around SynHydro's KirschGenerator.

This module provides a KirschBorgWrapper class that adapts SynHydro's Kirsch
nonparametric bootstrap generator to work with MOEA-FIND's optimization framework.

Two decision variable modes are supported:
  - "index": DVs [0,1] are mapped to discrete year indices via floor(p * n_years_hist),
    passed to KirschGenerator.generate_from_indices(indices).
  - "residual": DVs [0,1] are mapped to standardized residuals via empirical CDF
    of historical residuals (per calendar month), passed to KirschGenerator.generate_from_residuals(residuals).

Both modes generate synthetic monthly streamflow traces that preserve the statistical
structure (autocorrelation, seasonal patterns, cross-site correlation) of the
historical data via Cholesky decomposition and normal-score transforms.

References
----------
- Kirsch, B.R., Characklis, G.W., and Zeff, H.B. (2013). Evaluating the impact
  of alternative hydro-climate scenarios on transfer agreements. Journal of Water
  Resources Planning and Management, 139(4), 396-406.
- SynHydro documentation: synhydro/methods/generation/nonparametric/kirsch.py
"""

from typing import Literal, Optional

import numpy as np


class KirschBorgWrapper:
    """Borg-compatible wrapper for SynHydro's KirschGenerator.

    Maps Borg decision variables [0,1] to either year indices or residuals,
    then generates synthetic monthly flows using the full Kirsch pipeline
    (Cholesky decomposition, normal-score transform, seasonal adjustment).

    Attributes:
        kirsch_gen: Fitted SynHydroKirschGenerator instance.
        mode: Decision variable mapping mode ("index" or "residual").
        n_years_out: Number of output years to generate.
        n_years_hist: Number of historical years (for index mapping).
        sorted_residuals: Sorted standardized residuals by month (for residual mode).
    """

    def __init__(
        self,
        kirsch_gen,
        mode: Literal["index", "residual"] = "index",
        n_years_out: int = 30,
    ):
        """Initialize KirschBorgWrapper.

        Args:
            kirsch_gen: Fitted SynHydro KirschGenerator instance.
                Must have been preprocessed and fitted before passing here.
            mode: Decision variable mapping mode.
                - "index": DVs map to year indices via floor(p * n_years_hist).
                - "residual": DVs map to standardized residuals via empirical CDF.
            n_years_out: Number of years to generate per call to generate().

        Raises:
            ValueError: If mode is not "index" or "residual".
            AttributeError: If kirsch_gen is not fitted (missing Z_h, mean_month, std_month).
        """
        if mode not in ("index", "residual"):
            raise ValueError(f"mode must be 'index' or 'residual', got {mode!r}")

        self.kirsch_gen = kirsch_gen
        self.mode = mode
        self.n_years_out = n_years_out

        # Validate that the Kirsch generator is fitted
        if not hasattr(kirsch_gen, "Z_h"):
            raise AttributeError(
                "KirschGenerator must be fitted before wrapping. "
                "Call kirsch_gen.fit(Q_obs) first."
            )
        if not hasattr(kirsch_gen, "mean_month"):
            raise AttributeError(
                "KirschGenerator must have mean_month attribute. "
                "Ensure fit() was called."
            )

        self.n_years_hist = kirsch_gen.n_historic_years
        self.n_sites = kirsch_gen.n_sites

        # For residual mode, pre-compute sorted residuals per month for CDF mapping
        self.sorted_residuals = None
        if mode == "residual":
            self._prepare_residuals()

    def _prepare_residuals(self) -> None:
        """Pre-compute sorted standardized residuals per month for CDF mapping.

        Standardized residuals are computed as:
            Z = (log(Q) - mean_month) / std_month

        for each month across all historical years. These are sorted to enable
        empirical CDF inversion: DV p maps to the p-quantile of residuals.
        """
        # Z_h has shape (n_years, 12, n_sites)
        self.sorted_residuals = np.sort(self.kirsch_gen.Z_h, axis=0)
        # sorted_residuals[i, m, s] = i-th sorted residual for month m, site s

    def _map_indices(self, dvs: np.ndarray) -> np.ndarray:
        """Map decision variables [0,1] to year indices [0, n_years_hist-1].

        Args:
            dvs: Flat array of DVs in [0,1], length must equal n_dvs.

        Returns:
            Array of year indices with same shape as dvs.
        """
        # floor(p * n_years_hist) gives indices in [0, n_years_hist-1]
        indices = np.floor(dvs * self.n_years_hist).astype(int)
        # Clip to ensure no out-of-bounds
        return np.clip(indices, 0, self.n_years_hist - 1)

    def _map_residuals(self, dvs: np.ndarray) -> np.ndarray:
        """Map decision variables [0,1] to standardized residuals via empirical CDF.

        For each output month m and each site s, DV p maps to the p-quantile
        of sorted historical residuals for that month/site.

        Args:
            dvs: Flat array of DVs in [0,1], length must equal n_dvs.

        Returns:
            Array of standardized residuals with shape (n_years_out, 12, n_sites).
        """
        dvs_reshaped = dvs.reshape(self.n_years_out, 12)
        residuals = np.zeros((self.n_years_out, 12, self.n_sites))

        for year in range(self.n_years_out):
            for month in range(12):
                for site in range(self.n_sites):
                    p = dvs_reshaped[year, month]
                    # Map p to the p-quantile of sorted residuals for this month/site
                    idx = int(np.clip(p * self.n_years_hist, 0, self.n_years_hist - 1))
                    residuals[year, month, site] = (
                        self.sorted_residuals[idx, month, site]
                    )

        return residuals

    @property
    def n_dvs(self) -> int:
        """Total number of decision variables needed.

        For index mode: (n_years_out + 1) * 12, because KirschGenerator needs
        one extra year for Dec-Jan cross-year correlation handling.

        For residual mode: n_years_out * 12.

        Multisite uses the same number (per-month decisions), with each
        decision affecting all sites together via the Kirsch correlation structure.
        """
        if self.mode == "index":
            return (self.n_years_out + 1) * 12
        return self.n_years_out * 12

    def generate(self, dvs: np.ndarray) -> np.ndarray:
        """Generate synthetic monthly flows from decision variables.

        Based on the wrapper's mode, maps DVs to either year indices or residuals,
        then calls the appropriate KirschGenerator method. Returns a 2D array of
        synthetic monthly flows.

        Args:
            dvs: Flat array of decision variables in [0,1].
                Length must equal self.n_dvs (n_years_out * 12).

        Returns:
            Synthetic monthly flows with shape (n_years_out, 12) for single-site
            or (n_years_out * 12, n_sites) depending on KirschGenerator output format.
            Reshaping may be needed by the caller depending on use case.

        Raises:
            ValueError: If len(dvs) != self.n_dvs.
            AttributeError: If KirschGenerator lacks required methods
                (generate_from_indices or generate_from_residuals).
        """
        if len(dvs) != self.n_dvs:
            raise ValueError(
                f"Expected {self.n_dvs} decision variables, got {len(dvs)}"
            )

        if self.mode == "index":
            indices = self._map_indices(dvs)
            # Reshape to (n_years_out + 1, 12) for KirschGenerator
            # The extra year is needed for Dec-Jan cross-year handling
            indices_2d = indices.reshape(self.n_years_out + 1, 12)

            if not hasattr(self.kirsch_gen, "generate_from_indices"):
                raise AttributeError(
                    "KirschGenerator must have generate_from_indices() method. "
                    "Ensure SynHydro version supports this."
                )
            synthetic = self.kirsch_gen.generate_from_indices(
                indices_2d, n_years=self.n_years_out
            )

        elif self.mode == "residual":
            residuals = self._map_residuals(dvs)
            # residuals has shape (n_years_out, 12, n_sites)

            # Call KirschGenerator.generate_from_residuals()
            if not hasattr(self.kirsch_gen, "generate_from_residuals"):
                raise AttributeError(
                    "KirschGenerator must have generate_from_residuals() method. "
                    "Ensure SynHydro version supports this."
                )
            synthetic = self.kirsch_gen.generate_from_residuals(residuals)

        # Ensure output is numpy array (Kirsch may return DataFrame)
        if hasattr(synthetic, "values"):
            synthetic = synthetic.values

        # Reshape to 2D (n_years_out, 12) for single-site
        if synthetic.ndim == 1:
            synthetic = synthetic.reshape(self.n_years_out, 12)
        elif synthetic.ndim == 2 and synthetic.shape[0] == self.n_years_out * 12:
            if self.n_sites == 1:
                synthetic = synthetic.reshape(self.n_years_out, 12)

        # SynHydro outputs calendar-year order (Jan=col 0, Dec=col 11).
        # The rest of the pipeline uses water-year order (Oct=col 0, Sep=col 11)
        # because prepare_data() aligns historical flows to Oct-Sep water years.
        # Roll columns to convert: Jan(0)..Dec(11) → Oct(0)..Sep(11).
        if synthetic.ndim == 2 and synthetic.shape[1] == 12:
            synthetic = np.roll(synthetic, 3, axis=1)

        return synthetic
