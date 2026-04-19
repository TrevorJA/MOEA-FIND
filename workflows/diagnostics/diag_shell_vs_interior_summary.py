"""Dimension-sweep summary for the shell-versus-interior diagnostic.

Reads every `outputs/diag_shell_vs_interior/k{K}/results.json` file and
produces a single summary figure that compares MOEA-FIND to the three
reference samplers across the dimensions tested. The figure is
Figure SI-2f in the supporting information.

Run after all individual k runs have completed:
    python scripts/diag_shell_vs_interior_summary.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_results(outputs_root: Path) -> list:
    """Load every available k{N}/results.json under outputs/diag_shell_vs_interior/."""
    rows = []
    pattern = re.compile(r"^k(\d+)$")
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir():
            continue
        m = pattern.match(child.name)
        if not m:
            continue
        k = int(m.group(1))
        rj = child / "results.json"
        if not rj.exists():
            print(f"  skip k={k}: missing results.json")
            continue
        data = json.loads(rj.read_text())
        rows.append((k, data))
    rows.sort(key=lambda t: t[0])
    return rows


def _sampler_keys(row: dict) -> list[tuple[str, str]]:
    """Return [(json_key, display_label), ...] for the reference samplers.

    Handles both feasibility shapes: cube runs write keys like
    ``uniform_in_cube`` / ``lhs_in_cube`` / ``sobol_in_cube``; ball runs
    use the ``_in_ball`` suffix. This helper lets the summary figure
    work transparently across shapes.
    """
    shape = row.get("feasible_shape", "ball")
    ref_suffix = f"_in_{shape}"
    return [
        ("MOEA-FIND", "MOEA-FIND"),
        (f"uniform{ref_suffix}", "uniform"),
        (f"lhs{ref_suffix}", "LHS"),
        (f"sobol{ref_suffix}", "Sobol"),
    ]


def plot_summary(rows: list, out_path: Path) -> None:
    ks = [k for k, _ in rows]
    # Sampler keys can differ row-by-row if the sweep mixes shapes. For
    # labels + colors we drive off the first row, which in practice is
    # consistent across a sweep.
    sampler_pairs = _sampler_keys(rows[0][1])
    samplers = [s for s, _ in sampler_pairs]
    labels = [lbl for _, lbl in sampler_pairs]
    color_list = ["tab:blue", "tab:gray", "tab:orange", "tab:green"]
    colors = dict(zip(samplers, color_list))

    # Four panels: mean distance to D*, interior fraction, orthant occupancy,
    # grid cell occupancy. Each panel plots one line per sampler as a function
    # of k, with markers at each tested dimensionality.
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    def _get_mean(row, sampler):
        return row[sampler]["dist_from_D*"]["mean"]

    def _get_interior(row, sampler):
        return row[sampler]["interior_fraction"]

    def _get_orthant(row, sampler):
        occ = row[sampler].get("orthant_occupancy", None)
        if not occ:
            return float("nan")
        return occ.get("fraction", float("nan"))

    def _get_grid(row, sampler):
        occ = row[sampler].get("grid_occupancy", None)
        if not occ:
            return float("nan")
        return occ.get("fraction", float("nan"))

    getters = [
        ("Mean L1 distance from $D^*$", _get_mean),
        ("Interior mass fraction", _get_interior),
        (f"Orthant occupancy ($2^k$ orthants)", _get_orthant),
        ("Grid cell coverage fraction", _get_grid),
    ]

    for ax, (title, getter) in zip(axes, getters):
        for sampler, label in zip(samplers, labels):
            ys = [getter(r[1], sampler) for r in rows]
            ax.plot(ks, ys, marker="o", color=colors[sampler], label=label)
        ax.set_xlabel("dimensionality k")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(ks)
    axes[0].legend(fontsize=8, loc="best")

    shape = rows[0][1].get("feasible_shape", "ball")
    radius = rows[0][1].get("feasible_radius", rows[0][1].get("ball_radius", 2.5))
    fig.suptitle(
        f"Shell versus interior diagnostic, dimension sweep "
        f"(k-{shape} half-width {radius}, anti-ideal at the positive corner)"
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def write_table(rows: list, out_path: Path) -> None:
    """Write a plain Markdown summary table to accompany the figure."""
    lines = [
        "| k | shape | sampler    | n archive | mean L1 to D* | std | interior frac | orthant frac | grid frac |",
        "|---|-------|------------|-----------|---------------|-----|---------------|--------------|-----------|",
    ]
    for k, data in rows:
        shape = data.get("feasible_shape", "ball")
        for key, name in _sampler_keys(data):
            s = data[key]
            d = s["dist_from_D*"]
            orth = s.get("orthant_occupancy", {}).get("fraction", float("nan"))
            grid = s.get("grid_occupancy", {}).get("fraction", float("nan"))
            lines.append(
                f"| {k} | {shape} | {name:<10} | {s['n']} | "
                f"{d['mean']:.3f} | {d['std']:.3f} | "
                f"{s['interior_fraction']:.3f} | {orth:.3f} | {grid:.3f} |"
            )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--outputs",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "diag_shell_vs_interior",
    )
    p.add_argument(
        "--figure",
        type=Path,
        default=PROJECT_ROOT / "figures" / "figSI_shell_interior_sweep.pdf",
    )
    p.add_argument(
        "--table",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "diag_shell_vs_interior" / "sweep_table.md",
    )
    args = p.parse_args()

    rows = load_results(args.outputs)
    if not rows:
        print("no k{N} results found")
        sys.exit(1)
    print(f"loaded {len(rows)} k-runs: {[k for k, _ in rows]}")

    plot_summary(rows, args.figure)
    write_table(rows, args.table)
    print(f"figure: {args.figure}")
    print(f"table:  {args.table}")


if __name__ == "__main__":
    main()
