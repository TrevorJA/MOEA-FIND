"""Project historical T-year blocks onto each registered metric set.

For motivation work (DD-04 review): given the historical Cannonsville
record, slide a T-year window across it and compute drought
characteristics on each block. Project the resulting cloud onto each
metric set in :data:`src.drought_metrics.PRESETS` and produce a 3D
scatter figure per preset, plus a per-metric range/coefficient-of-variation
table.

The figures answer the question "which drought metrics reflect truly
uncertain properties of the historical record?" — a metric whose
historical-block cloud is tightly clustered is essentially constant
across the basin's recorded history and therefore makes a weak axis
for drought-hazard exploration. A metric whose cloud spreads broadly
is one MOEA-FIND can usefully explore.

Reproducible: deterministic (sliding-block extraction is deterministic),
seed not used, configuration logged to ``config.json``. Output:
``outputs/02_calibration/metric_blocks/``. Plotting lives in
``workflows/02_calibration/plots/metric_blocks.py``.

Usage::

    python workflows/02_calibration/metric_blocks.py \\
        --T-years 20 --stride 1 --ssi 3
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.drought_metrics import (  # noqa: E402
    PRESETS,
    REGISTRY,
    metric_labels,
    metric_names,
    resolve_metric_set,
)
from src.experiment_utils import prepare_data, compute_historical_ssi_chars  # noqa: E402
from src.historical_blocks import resample_historical_blocks  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402
from src.objectives import (  # noqa: E402
    compute_ssi,
    compute_ssi_drought_characteristics,
    flows_to_series,
)

STAGE = "02_calibration"
DRIVER = "metric_blocks"


def compute_block_chars(
    monthly_1d: np.ndarray,
    T_years: int,
    stride: int,
    ssi_calc,
) -> pd.DataFrame:
    """Return a DataFrame with one row per T-year block and one column
    per chars-dict key produced by
    :func:`compute_ssi_drought_characteristics`.

    SSI is computed by transforming each block through the prefitted
    ``ssi_calc`` so that all blocks share the same calibrated
    distribution (this is the convention the manuscript uses; see
    DD-11 lock-in note in :func:`src.historical_blocks.compute_historical_block_chars`).
    The flows array is passed through so the chars dict includes
    ``q10_flow_neg`` for trace-level metrics.
    """
    blocks = resample_historical_blocks(monthly_1d, T_years, stride)
    rows: List[Dict[str, float]] = []
    for i, blk in enumerate(blocks):
        ssi = ssi_calc.transform(flows_to_series(blk, start_date="2100-01-01"))
        chars = compute_ssi_drought_characteristics(ssi, monthly_flows=blk)
        # Drop n_events (integer count) so the DataFrame is all numeric.
        chars = {k: v for k, v in chars.items() if k != "n_events"}
        rows.append(chars)
    return pd.DataFrame(rows)


def plot_metric_set_3d(
    chars_df: pd.DataFrame,
    preset_name: str,
    metric_set,
    fig_path: Path,
    *,
    historical_aggregate: dict | None = None,
) -> None:
    """One 3D scatter plot per metric set; one point per T-block."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers projection)

    if len(metric_set) < 3:
        return  # all current presets are 3D; skip if a future preset is shorter

    pts = np.array([[m.extract(row) for m in metric_set] for row in chars_df.to_dict("records")])

    fig = plt.figure(figsize=(7.0, 6.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        pts[:, 0], pts[:, 1], pts[:, 2],
        c=np.arange(len(pts)),
        cmap="viridis",
        s=24, alpha=0.85, edgecolors="black", linewidths=0.3,
    )
    if historical_aggregate is not None:
        ha = np.array([m.extract(historical_aggregate) for m in metric_set])
        ax.scatter(
            [ha[0]], [ha[1]], [ha[2]],
            marker="*", s=220, c="red", edgecolors="black", linewidths=0.6,
            label="full historical record",
        )
        ax.legend(loc="upper left", fontsize=8)

    m0, m1, m2 = metric_set[0], metric_set[1], metric_set[2]
    ax.set_xlabel(f"{m0.label}\n({m0.units})", fontsize=9)
    ax.set_ylabel(f"{m1.label}\n({m1.units})", fontsize=9)
    ax.set_zlabel(f"{m2.label}\n({m2.units})", fontsize=9)
    ax.set_title(
        f"Preset: {preset_name}\n"
        f"{len(pts)} historical T-year blocks (colour = chronological order)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(fig_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(fig_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def per_metric_summary(chars_df: pd.DataFrame) -> pd.DataFrame:
    """Range, std, and coefficient of variation for every registry metric.

    A high CV means the metric varies meaningfully across historical
    T-year blocks. A low CV means the metric is roughly constant across
    the historical record and is therefore a weak axis for drought-
    hazard exploration.
    """
    rows = []
    for name, m in REGISTRY.items():
        col_name = "q10_flow_neg" if name == "q10_flow" else name
        if col_name not in chars_df.columns:
            continue
        vals = chars_df[col_name].astype(float).values
        if vals.size == 0:
            continue
        mean = float(np.mean(vals))
        std = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
        cv = float(std / abs(mean)) if abs(mean) > 1e-12 else float("nan")
        rows.append({
            "metric": name,
            "label": m.label,
            "units": m.units,
            "mean": mean,
            "std": std,
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
            "range": float(np.max(vals) - np.min(vals)),
            "cv": cv,
        })
    return pd.DataFrame(rows).sort_values("cv", ascending=False).reset_index(drop=True)


def per_preset_summary(
    chars_df: pd.DataFrame,
    presets: Dict[str, tuple],
) -> pd.DataFrame:
    """For each preset, the combined volume and per-axis ranges of the
    historical-block cloud (a rough proxy for "how much hazard space
    can a structured sampler hope to populate from the historical
    envelope alone")."""
    rows = []
    for preset_name in presets:
        ms = resolve_metric_set(preset_name)
        ranges = []
        for m in ms:
            col_name = "q10_flow_neg" if m.name == "q10_flow" else m.name
            if col_name not in chars_df.columns:
                ranges.append(float("nan"))
                continue
            v = chars_df[col_name].astype(float).values
            ranges.append(float(np.max(v) - np.min(v)) if v.size else 0.0)
        rows.append({
            "preset": preset_name,
            "metrics": ", ".join(metric_names(ms)),
            "axis_ranges": ", ".join(f"{r:.3g}" for r in ranges),
            "log_volume_proxy": float(
                np.sum(np.log(np.maximum(np.array(ranges), 1e-12)))
            ),
        })
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--T-years", type=int, default=20,
                   help="Length of each T-year block (default 20).")
    p.add_argument("--stride", type=int, default=1,
                   help="Block stride in years (default 1 → maximum overlap).")
    p.add_argument("--ssi", type=int, default=3,
                   choices=[1, 3, 6, 12],
                   help="SSI accumulation period in months.")
    p.add_argument("--cache-dir", type=Path,
                   default=PROJECT_ROOT / "outputs" / "data_cache")
    args = p.parse_args()

    out_dir = stage_output_dir(STAGE, DRIVER)
    print(f"[02/metric_blocks] output_dir={out_dir}")

    # --- Data ---
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(args.cache_dir)
    n_hist_years = monthly_2d.shape[0]
    n_blocks = (n_hist_years - args.T_years) // args.stride + 1
    print(f"[diag] historical: {n_hist_years} water years, "
          f"{n_blocks} blocks of T={args.T_years} years (stride={args.stride})")

    # --- Prefit SSI calculator on full historical record ---
    _, ssi_calc, full_hist_chars = compute_historical_ssi_chars(monthly_1d, args.ssi)
    # Augment full-record chars with trace-level extras for the historical star.
    ssi_full_series, _ = compute_ssi(monthly_1d, timescale=args.ssi)
    full_hist_chars = compute_ssi_drought_characteristics(
        ssi_full_series, monthly_flows=monthly_1d
    )

    # --- Per-block chars (one row per block, all keys present) ---
    print(f"[diag] computing per-block drought characteristics ...")
    chars_df = compute_block_chars(
        monthly_1d, args.T_years, args.stride, ssi_calc
    )
    chars_df.to_csv(out_dir / "block_chars.csv", index=False)
    print(f"[diag] wrote {out_dir / 'block_chars.csv'} ({len(chars_df)} rows, "
          f"{len(chars_df.columns)} columns)")

    # --- Per-metric summary table ---
    metric_summary = per_metric_summary(chars_df)
    metric_summary.to_csv(out_dir / "per_metric_summary.csv", index=False)
    print(f"[diag] wrote {out_dir / 'per_metric_summary.csv'}")
    print(metric_summary.to_string(index=False))

    # --- Per-preset summary table ---
    preset_summary = per_preset_summary(chars_df, PRESETS)
    preset_summary.to_csv(out_dir / "per_preset_summary.csv", index=False)
    print(f"\n[diag] wrote {out_dir / 'per_preset_summary.csv'}")
    print(preset_summary.to_string(index=False))

    # Dump full historical chars for the plotting driver to render the
    # red-star reference point.
    import pickle
    (out_dir / "full_hist_chars.pkl").write_bytes(pickle.dumps(full_hist_chars))

    # --- Reproducibility: dump config ---
    cfg = {
        "script": "diag_metric_set_historical_blocks.py",
        "T_years": args.T_years,
        "stride": args.stride,
        "ssi": args.ssi,
        "n_hist_years": int(n_hist_years),
        "n_blocks": int(len(chars_df)),
        "presets": {k: list(v) for k, v in PRESETS.items()},
        "registry_keys": sorted(REGISTRY.keys()),
    }
    (out_dir / "config.json").write_text(json.dumps(cfg, indent=2))
    print(f"\n[diag] wrote {out_dir / 'config.json'}")
    print(f"[diag] done.")


if __name__ == "__main__":
    main()
