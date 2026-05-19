"""Per-mode quick-look plot for the wrapper-mode ablation.

Renders a single drought-space scatter for one ``(mode, slug)`` pair under
``figures/04_moea_find_single_site/wrapper_mode_ablation/<mode>/<slug>/``.
Use the companion ``plots/wrapper_mode_compare.py`` for the SI panel set.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402

STAGE = "04_moea_find_single_site"
DRIVER = "wrapper_mode_ablation"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--mode", required=True, choices=["index", "residual"])
    p.add_argument("--slug", required=True)
    args = p.parse_args()

    in_dir = stage_output_dir(STAGE, DRIVER, f"{args.mode}/{args.slug}", create=False)
    results_path = in_dir / "results.json"
    if not results_path.exists():
        sys.exit(f"missing {results_path} -- run the compute driver first")
    r = json.loads(results_path.read_text())
    dm = np.array(r.get("drought_metrics", []), dtype=float)
    if dm.size == 0:
        print("  no Pareto solutions; nothing to render.")
        return

    fig_dir = stage_figure_dir(STAGE, DRIVER, f"{args.mode}/{args.slug}")
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    ax.scatter(dm[:, 0], dm[:, 1], s=12, alpha=0.85, color="#c05621",
               label=f"Pareto (n={len(dm)})")
    ax.set_xlabel("Mean duration (months)")
    ax.set_ylabel("Mean avg. severity")
    ax.set_title(f"Wrapper-mode ablation: mode={args.mode}\n{args.slug}")
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()
    out = fig_dir / "drought_space.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots/04/wrapper_mode_ablation] wrote {out}")


if __name__ == "__main__":
    main()
