"""Multi-panel 2D + 3D scatter of a first-event MMBorg archive vs history.

The stock stage-04 plots driver only renders a single 2-objective 2D
scatter + one 3D triplet, and has NO historical overlay for first-event
presets (stage 04 skips historical_block_chars.npz for them). This
diagnostic fills that gap for the windowed/first-event coarse archive:

  * Full pairwise 2D scatter MATRIX over all K objectives (lower
    triangle = scatter, diagonal = marginal histograms), synthetic
    Pareto cloud vs a historical first-event cloud, with the anti-ideal
    D* reference per axis.
  * Several 3D scatter views over informative objective triplets,
    same synthetic-vs-historical overlay + D* marker.

Historical reference (consistent with the first-event framing): every
overlapping T-year window of the historical record is reduced to the
characteristics of ITS first critical SSI-3 drought (same extractor as
the stage-04 objectives, incl. the 4.5 severity clip), giving a
historical cloud directly comparable to the synthetic archive. The
single full-record first-event point is also marked.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.experiment import prepare_data  # noqa: E402
from src.metrics.objectives import (  # noqa: E402
    make_ssi_calculator, flows_to_series, compute_ssi_drought_characteristics,
)
from src.metrics.drought_metrics import REGISTRY  # noqa: E402
from src.hydrology.historical_blocks import resample_historical_blocks  # noqa: E402

SYN_KW = dict(s=10, alpha=0.35, c="#1f77b4", edgecolors="none", label="synthetic Pareto")
HIST_KW = dict(s=26, alpha=0.85, c="#d62728", marker="x", label="historical T-blocks")
FULL_KW = dict(s=140, c="#000000", marker="*", label="historical full record")


def _historical_first_event_cloud(objective_keys, ssi_timescale, n_years):
    """Per-block historical first-event metric rows (n_blocks_with_event x K)."""
    cache = PROJECT_ROOT / "outputs" / "data_cache"
    _, monthly_1d = prepare_data(cache)
    ssi = make_ssi_calculator(timescale=ssi_timescale)
    ssi.fit(flows_to_series(monthly_1d, start_date="1950-10-01"))
    extracts = [REGISTRY[k].extract for k in objective_keys]

    def _row(flat):
        s = ssi.transform(flows_to_series(np.asarray(flat, float),
                                          start_date="2100-01-01"))
        ch = compute_ssi_drought_characteristics(s, monthly_flows=flat)
        if int(ch.get("first_event_present", 0)) == 0:
            return None
        return [float(e(ch)) for e in extracts]

    rows = []
    for blk in resample_historical_blocks(monthly_1d, T_years=n_years, stride=1):
        r = _row(blk)
        if r is not None:
            rows.append(r)
    full = _row(monthly_1d)
    return (np.asarray(rows, float) if rows else np.empty((0, len(objective_keys))),
            np.asarray(full, float) if full is not None else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--ssi", type=int, default=3)
    args = ap.parse_args()

    R = json.loads(Path(args.results).read_text())
    syn = np.asarray(R["drought_metrics"], float)          # (N, K)
    keys = list(R["objective_keys"])
    K = len(keys)
    syn = syn[:, :K]
    dstar = np.asarray(R.get("anti_ideal", [np.nan] * K), float)
    n_years = int(R.get("n_years_out", 10))
    labels = [f"{REGISTRY[k].label}\n[{REGISTRY[k].units}]" for k in keys]
    short = [REGISTRY[k].label for k in keys]

    print(f"[fe-scatter] synthetic archive: {syn.shape[0]} x {K}", flush=True)
    hist, full = _historical_first_event_cloud(keys, args.ssi, n_years)
    print(f"[fe-scatter] historical first-event cloud: {hist.shape[0]} blocks; "
          f"full-record point={'yes' if full is not None else 'no'}", flush=True)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- A) pairwise 2D scatter matrix --------------------------------
    fig, axes = plt.subplots(K, K, figsize=(3.1 * K, 3.0 * K))
    for i in range(K):
        for j in range(K):
            ax = axes[i, j]
            if i == j:
                ax.hist(syn[:, i], bins=30, color=SYN_KW["c"], alpha=0.55,
                        density=True)
                if hist.shape[0]:
                    ax.hist(hist[:, i], bins=20, color=HIST_KW["c"],
                            alpha=0.55, density=True)
                if np.isfinite(dstar[i]):
                    ax.axvline(dstar[i], color="k", ls="--", lw=1)
                ax.set_yticks([])
            elif i > j:
                if hist.shape[0]:
                    ax.scatter(hist[:, j], hist[:, i], **HIST_KW)
                ax.scatter(syn[:, j], syn[:, i], **SYN_KW)
                if full is not None:
                    ax.scatter([full[j]], [full[i]], **FULL_KW)
                if np.isfinite(dstar[j]):
                    ax.axvline(dstar[j], color="k", ls="--", lw=0.8, alpha=0.6)
                if np.isfinite(dstar[i]):
                    ax.axhline(dstar[i], color="k", ls="--", lw=0.8, alpha=0.6)
            else:
                ax.axis("off")
            if j == 0 and i != 0:
                ax.set_ylabel(labels[i], fontsize=8)
            if i == K - 1:
                ax.set_xlabel(labels[j], fontsize=8)
            ax.tick_params(labelsize=7)
    # Explicit Line2D legend proxies (do NOT splat the scatter kwargs:
    # edgecolors/c/s are Collection kwargs Line2D rejects).
    h = [
        plt.Line2D([], [], linestyle="", marker="o", markersize=8,
                   color=SYN_KW["c"]),
        plt.Line2D([], [], linestyle="", marker="x", markersize=8,
                   color=HIST_KW["c"]),
        plt.Line2D([], [], linestyle="", marker="*", markersize=12,
                   color=FULL_KW["c"]),
    ]
    fig.legend(h, ["synthetic Pareto", "historical T-blocks",
                   "historical full record"],
               loc="upper center", ncol=3, fontsize=10)
    fig.suptitle(f"First-event MMBorg archive vs historical "
                 f"(K={K}, N_syn={syn.shape[0]}, dashed = anti-ideal D*)",
                 y=0.995, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p2d = out / "scatter_matrix_2d.pdf"
    fig.savefig(p2d, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fe-scatter] wrote {p2d}", flush=True)

    # ---- B) 3D scatter views over informative triplets ----------------
    triplets = list(itertools.combinations(range(K), 3))
    # Cap at 6 panels; prefer triplets that include severity (axis 0).
    triplets.sort(key=lambda t: (0 not in t, t))
    triplets = triplets[:6]
    ncol = 3
    nrow = int(np.ceil(len(triplets) / ncol))
    fig = plt.figure(figsize=(6.2 * ncol, 5.4 * nrow))
    for idx, (a, b, c) in enumerate(triplets):
        ax = fig.add_subplot(nrow, ncol, idx + 1, projection="3d")
        if hist.shape[0]:
            ax.scatter(hist[:, a], hist[:, b], hist[:, c],
                       s=22, c=HIST_KW["c"], marker="x", alpha=0.85)
        ax.scatter(syn[:, a], syn[:, b], syn[:, c],
                   s=8, c=SYN_KW["c"], alpha=0.30, edgecolors="none")
        if full is not None:
            ax.scatter([full[a]], [full[b]], [full[c]], s=120, c="k",
                       marker="*")
        if np.all(np.isfinite(dstar[[a, b, c]])):
            ax.scatter([dstar[a]], [dstar[b]], [dstar[c]], s=90, c="k",
                       marker="X")
        ax.set_xlabel(short[a], fontsize=7)
        ax.set_ylabel(short[b], fontsize=7)
        ax.set_zlabel(short[c], fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_title(f"{short[a]} / {short[b]} / {short[c]}", fontsize=8)
    fig.suptitle("First-event MMBorg archive vs historical — 3D views "
                 "(blue=synthetic, red x=historical blocks, "
                 "★=full-record, X=anti-ideal D*)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p3d = out / "scatter_3d_views.pdf"
    fig.savefig(p3d, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fe-scatter] wrote {p3d}", flush=True)
    print(f"[fe-scatter] DONE -> {out}", flush=True)


if __name__ == "__main__":
    main()
