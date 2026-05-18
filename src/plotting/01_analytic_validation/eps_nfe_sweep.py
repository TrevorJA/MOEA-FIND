"""Render the eps x NFE heatmap from the cached aggregate.json."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from src.io_paths.paths import stage_figure_dir, stage_output_dir  # noqa: E402
from src.plotting.analytic import fig3_eps_nfe_heatmap  # noqa: E402

STAGE = "01_analytic_validation"
DRIVER = "eps_nfe_sweep"


def main():
    out_dir = stage_output_dir(STAGE, DRIVER, create=False)
    agg_path = out_dir / "aggregate.json"
    if not agg_path.exists():
        sys.exit(f"missing {agg_path} — run the compute driver in --mode aggregate first")
    agg = json.loads(agg_path.read_text())["aggregated"]

    fig_dir = stage_figure_dir(STAGE, DRIVER)
    out_pdf = fig_dir / "eps_nfe_sweep.pdf"
    fig = fig3_eps_nfe_heatmap(agg)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"[plots/01/eps_nfe_sweep] {out_pdf}")


if __name__ == "__main__":
    main()
