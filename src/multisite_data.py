"""Multi-site data loading and generator fitting for DRB policy re-evaluation.

Loads multi-site historical gage flows and catchment inflows from Pywr-DRB
datasets, identifies which nodes are Kirsch-generated vs KDE-regressed,
and fits SynHydro's KirschGenerator and NowakDisaggregator on the
appropriate node subsets.

Follows the same data loading and fitting patterns used by
StochasticExploratoryExperiment/methods/generate.py, but self-contained
within MOEA-FIND (no edits to external repos).
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple, List

import numpy as np
import pandas as pd

from synhydro.methods.generation.nonparametric.kirsch import KirschGenerator
from synhydro.methods.disaggregation.temporal.nowak import NowakDisaggregator

from pywrdrb.path_manager import get_pn_object
from pywrdrb.pywr_drb_node_data import (
    immediate_downstream_nodes_dict,
    downstream_node_lags,
)

# Default baseline dataset (same as StochasticExploratoryExperiment)
DEFAULT_BASELINE = "pub_nhmv10_BC_withObsScaled"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_pywrdrb_gage_flow(
    baseline_dataset: str = DEFAULT_BASELINE,
    period: str = "full",
) -> pd.DataFrame:
    """Load multi-site daily gage flows from Pywr-DRB data directory.

    Args:
        baseline_dataset: Name of the flow dataset directory under
            ``pywrdrb/data/flows/``.
        period: ``"full"`` for the entire record, or ``"baseline"``
            for 1980-2019.

    Returns:
        DataFrame with daily datetime index, site name columns, units MGD.
    """
    pn = get_pn_object()
    fpath = Path(pn.sc.get(f"flows/{baseline_dataset}")) / "gage_flow_mgd.csv"
    Q = pd.read_csv(fpath, index_col=0, parse_dates=True)
    Q.index = pd.to_datetime(Q.index)

    if period == "baseline":
        Q = Q.loc["1980-01-01":"2019-12-31"]

    # Replace zeros with NaN (physically unrealistic for gage flows)
    Q.replace(0, np.nan, inplace=True)
    return Q


def load_pywrdrb_catchment_inflow(
    baseline_dataset: str = DEFAULT_BASELINE,
    period: str = "full",
) -> pd.DataFrame:
    """Load multi-site daily catchment inflows from Pywr-DRB data directory.

    Args:
        baseline_dataset: Name of the flow dataset directory.
        period: ``"full"`` or ``"baseline"`` (1980-2019).

    Returns:
        DataFrame with daily datetime index, site name columns, units MGD.
    """
    pn = get_pn_object()
    fpath = Path(pn.sc.get(f"flows/{baseline_dataset}")) / "catchment_inflow_mgd.csv"
    Q = pd.read_csv(fpath, index_col=0, parse_dates=True)
    Q.index = pd.to_datetime(Q.index)

    if period == "baseline":
        Q = Q.loc["1980-01-01":"2019-12-31"]

    return Q


# ---------------------------------------------------------------------------
# Node classification
# ---------------------------------------------------------------------------

def get_kirsch_sites(Q_gage: pd.DataFrame) -> List[str]:
    """Return sites generated directly by Kirsch (non-USGS-gage nodes).

    Filters out USGS gage ID nodes (names starting with ``'0'``) and
    ``'delTrenton'`` (zeroed by convention). Follows the
    ``pywrdrb_nodes_to_generate`` pattern from
    StochasticExploratoryExperiment.

    Returns:
        List of ~22 major DRB node names present in *Q_gage*.
    """
    pywrdrb_nodes = list(immediate_downstream_nodes_dict.keys())
    kirsch_nodes = [n for n in pywrdrb_nodes if not n.startswith("0")]
    if "delTrenton" in kirsch_nodes:
        kirsch_nodes.remove("delTrenton")

    # Only keep nodes that exist in the gage flow data
    available = set(Q_gage.columns)
    return [n for n in kirsch_nodes if n in available]


def get_kde_regression_sites(Q_gage: pd.DataFrame) -> List[str]:
    """Return sites that need KDE regression (USGS gage ID nodes).

    These are nodes whose names start with ``'0'`` in the
    ``immediate_downstream_nodes_dict``.

    Returns:
        List of ~8 USGS gage node names present in *Q_gage*.
    """
    pywrdrb_nodes = list(immediate_downstream_nodes_dict.keys())
    kde_nodes = [n for n in pywrdrb_nodes if n.startswith("0")]

    available = set(Q_gage.columns)
    return [n for n in kde_nodes if n in available]


def get_kde_pairs(
    kirsch_sites: List[str],
    kde_sites: List[str],
) -> List[Tuple[str, str]]:
    """Return (upstream, downstream) pairs for KDE regression.

    For each Kirsch-generated upstream node whose immediate downstream
    is in *kde_sites*, yield the pair.

    Returns:
        List of ``(upstream_name, downstream_name)`` tuples.
    """
    pairs = []
    kde_set = set(kde_sites)
    for upstream in kirsch_sites:
        downstream = immediate_downstream_nodes_dict.get(upstream)
        if downstream in kde_set:
            pairs.append((upstream, downstream))
    return pairs


# ---------------------------------------------------------------------------
# Generator fitting
# ---------------------------------------------------------------------------

def fit_multisite_generators(
    Q_gage: pd.DataFrame,
    kirsch_sites: List[str],
) -> Tuple[KirschGenerator, NowakDisaggregator]:
    """Fit KirschGenerator (monthly) and NowakDisaggregator (daily) on
    multi-site gage flows.

    Steps:
        1. Filter *Q_gage* to *kirsch_sites* columns.
        2. Fit ``KirschGenerator(generate_using_log_flow=True)`` via
           ``.preprocessing()`` + ``.fit()`` on the filtered DataFrame.
        3. Fit ``NowakDisaggregator`` on the same filtered daily DataFrame.

    Args:
        Q_gage: Daily gage flow DataFrame (all sites, MGD).
        kirsch_sites: Column names to include in generator fitting.

    Returns:
        ``(kirsch_gen, nowak_disagg)`` — fitted generator and disaggregator.
    """
    Q_kirsch = Q_gage[kirsch_sites].copy()

    # Drop rows with any NaN to get a clean fitting window
    Q_kirsch = Q_kirsch.dropna()

    print(f"[multisite_data] Fitting Kirsch on {len(kirsch_sites)} sites, "
          f"{len(Q_kirsch)} days ({Q_kirsch.index[0].date()} to "
          f"{Q_kirsch.index[-1].date()})")

    kirsch_gen = KirschGenerator(generate_using_log_flow=True)
    kirsch_gen.preprocessing(Q_kirsch)
    kirsch_gen.fit()

    print(f"[multisite_data] Kirsch fitted: {kirsch_gen.n_sites} sites, "
          f"{kirsch_gen.n_historic_years} years")

    nowak_disagg = NowakDisaggregator()
    nowak_disagg.preprocessing(Q_kirsch)
    nowak_disagg.fit()

    print(f"[multisite_data] Nowak disaggregator fitted")

    return kirsch_gen, nowak_disagg
