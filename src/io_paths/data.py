"""Data loading utilities for MOEA-FIND."""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def load_usgs_daily(
    site_id: str = "01423000",
    start_date: str = "1950-10-01",
    end_date: str = "2023-09-30",
    cache_dir: Optional[Path] = None,
) -> pd.Series:
    """Load daily mean streamflow from USGS NWIS.

    Args:
        site_id: USGS gauge ID. Default is West Branch Delaware at Walton
            (Cannonsville inflow).
        start_date: Start date (YYYY-MM-DD).
        end_date: End date (YYYY-MM-DD).
        cache_dir: If provided, cache the downloaded data as CSV.

    Returns:
        pd.Series with DatetimeIndex and daily flow in cfs.
    """
    if cache_dir is not None:
        cache_path = Path(cache_dir) / f"usgs_{site_id}_daily.csv"
        if cache_path.exists():
            df = pd.read_csv(cache_path, index_col=0, parse_dates=True)
            return df["flow_cfs"]

    import urllib.request
    url = (
        f"https://waterservices.usgs.gov/nwis/dv/"
        f"?format=rdb&sites={site_id}"
        f"&startDT={start_date}&endDT={end_date}"
        f"&parameterCd=00060&statCd=00003"
    )
    resp = urllib.request.urlopen(url, timeout=30)
    lines = resp.read().decode().split("\n")

    # Parse RDB format: skip comments (#) and format line (5s...)
    data_lines = [
        l for l in lines
        if l and not l.startswith("#") and not l.startswith("5s")
    ]
    header = data_lines[0].split("\t")
    rows = [l.split("\t") for l in data_lines[1:] if l.strip()]

    dates = []
    flows = []
    for row in rows:
        if len(row) >= 4:
            try:
                dates.append(pd.Timestamp(row[2]))
                flows.append(float(row[3]))
            except (ValueError, IndexError):
                continue

    series = pd.Series(flows, index=pd.DatetimeIndex(dates), name="flow_cfs")
    series = series.sort_index()

    if cache_dir is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        series.to_frame().to_csv(cache_path)

    return series


def daily_to_monthly(daily: pd.Series) -> pd.Series:
    """Aggregate daily streamflow to monthly means.

    Args:
        daily: Daily flow series.

    Returns:
        Monthly mean flow series with period index.
    """
    monthly = daily.resample("MS").mean()
    # Drop months with insufficient data (< 25 days)
    counts = daily.resample("MS").count()
    monthly = monthly[counts >= 25]
    return monthly
