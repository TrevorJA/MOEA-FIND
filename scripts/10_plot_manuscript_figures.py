"""Script 10 — Manuscript figure assembler.

Reads the JSON/NPZ outputs of all preceding scripts (01-09) and writes the
final manuscript + SI figure set to ``figures/``. Each generated figure is
pinned to a specific section of the manuscript draft. The script is
deliberately *idempotent and stateless*: it regenerates any figure whose
upstream output exists and skips (with a warning) any figure whose upstream
data is missing. All plotting is delegated to src.plotting.

Outputs
-------
figures/fig01_manhattan_concept.pdf   §3.1   (no upstream data required)
figures/fig02_2d_analytic.pdf         §5.1   ← outputs/exp01_analytic_2d/pareto.npz
figures/fig02_3d_analytic.pdf         §5.2   ← outputs/exp02_analytic_3d/pareto.npz
figures/fig03_eps_nfe_sweep.pdf       §5.3   ← outputs/exp03_eps_nfe_sweep/aggregate.json
figures/fig05_kirsch_pareto.pdf       §6.1   ← outputs/exp04_kirsch_single_site/<mode>/results.json
figures/fig06_plausibility.pdf        §6.2   ← outputs/exp04_kirsch_single_site/<mode>/results.json
figures/fig07_coverage_comparison.pdf §6.3   ← outputs/exp05_kirsch_library/, outputs/exp06_library_subsample/
figures/fig08_drb_multisite.pdf       §7.1   ← outputs/exp08_drb_multisite/ (pending)
figures/fig09_policy_reeval.pdf       §7.3   ← outputs/exp09_drb_policy_reeval/ (pending)
figures/figSI01_hyperplane.pdf        SI-1
figures/figSI04_convergence.pdf       SI-4

Run:
    python scripts/10_plot_manuscript_figures.py          # produce everything available
    python scripts/10_plot_manuscript_figures.py --only fig03 fig07   # subset
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Callable, Dict, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Analytic drivers run without SynHydro installed.
from unittest.mock import MagicMock  # noqa: E402
_stub = MagicMock()
sys.modules.setdefault("synhydro", _stub)
sys.modules.setdefault("synhydro.droughts", _stub.droughts)
sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

import matplotlib.pyplot as plt  # noqa: E402

from src.plotting.analytic import (  # noqa: E402
    fig1_manhattan_concept,
    fig2_2d_coverage_comparison,
    fig2_3d_projections,
    fig3_eps_nfe_heatmap,
    fig_si_hyperplane_check,
)
from src.plotting.coverage import fig7_coverage_comparison  # noqa: E402
from src.plotting.style import COLORS  # noqa: E402


OUTPUTS = PROJECT_ROOT / "outputs"
FIGURES = PROJECT_ROOT / "figures"


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  wrote {path.relative_to(PROJECT_ROOT)}")


def _skip(name: str, reason: str) -> None:
    print(f"  SKIP {name}: {reason}")


# -----------------------------------------------------------------------------
# Figure factories — each one returns True on success, False if upstream missing.
# -----------------------------------------------------------------------------
def fig01_concept() -> bool:
    """§3.1 — Manhattan-norm conceptual figure (no upstream data)."""
    fig = fig1_manhattan_concept()
    _save(fig, FIGURES / "fig01_manhattan_concept.pdf")
    return True


def fig02_2d_analytic() -> bool:
    """§5.1 — 2D analytic Pareto vs QMC baselines."""
    npz = OUTPUTS / "exp01_analytic_2d" / "pareto.npz"
    if not npz.exists():
        _skip("fig02_2d", f"missing {npz.relative_to(PROJECT_ROOT)}")
        return False
    data = np.load(npz)
    fig = fig2_2d_coverage_comparison(data["dvs"], seed=42)
    _save(fig, FIGURES / "fig02_2d_analytic.pdf")
    return True


def fig02_3d_analytic() -> bool:
    """§5.2 — 3D analytic Pareto and simplex projections."""
    npz = OUTPUTS / "exp02_analytic_3d" / "pareto.npz"
    if not npz.exists():
        _skip("fig02_3d", f"missing {npz.relative_to(PROJECT_ROOT)}")
        return False
    data = np.load(npz)
    fig = fig2_3d_projections(data["dvs"])
    _save(fig, FIGURES / "fig02_3d_analytic.pdf")
    return True


def fig03_eps_nfe() -> bool:
    """§5.3 — epsilon × NFE sensitivity heatmap."""
    j = OUTPUTS / "exp03_eps_nfe_sweep" / "aggregate.json"
    if not j.exists():
        _skip("fig03", f"missing {j.relative_to(PROJECT_ROOT)}")
        return False
    agg = json.loads(j.read_text())["aggregated"]
    fig = fig3_eps_nfe_heatmap(agg)
    _save(fig, FIGURES / "fig03_eps_nfe_sweep.pdf")
    return True


def fig05_kirsch_pareto() -> bool:
    """§6.1 — single-site Kirsch Pareto in drought space."""
    from src.plotting.drought_space import plot_scatter_with_marginals
    for mode in ("residual", "index"):
        j = OUTPUTS / "exp04_kirsch_single_site" / mode / "results.json"
        if not j.exists():
            continue
        result = json.loads(j.read_text())
        pareto = np.asarray(result.get("pareto_drought_chars") or [])
        if pareto.size == 0:
            continue
        fig = plot_scatter_with_marginals(
            pareto, title=f"§6.1 Kirsch Pareto ({mode})",
            objective_labels=("Mean duration (months)", "Mean avg. severity"),
        )
        _save(fig, FIGURES / f"fig05_kirsch_pareto_{mode}.pdf")
        return True
    _skip("fig05", "no script 04 results found")
    return False


def fig06_plausibility() -> bool:
    """§6.2 — plausibility diagnostics (acf, FDC, seasonal cycle)."""
    # Script 04 --plot already writes individual plausibility PDFs.
    # Here we only stitch them into a multi-panel figure if the raw traces
    # are stored alongside results.json. Placeholder until the full trace
    # archive format is finalized.
    _skip("fig06", "full plausibility panel assembled inside script 04 --plot")
    return False


def fig07_coverage() -> bool:
    """§6.3 — headline coverage comparison (Fig 7)."""
    moea_json = OUTPUTS / "exp04_kirsch_single_site" / "residual" / "results.json"
    lib_json = OUTPUTS / "exp05_kirsch_library" / "characteristics.json"
    lhs_json = OUTPUTS / "exp06_library_subsample" / "subsample_lhs.json"
    sobol_json = OUTPUTS / "exp06_library_subsample" / "subsample_sobol.json"
    if not (moea_json.exists() and lib_json.exists() and lhs_json.exists()):
        _skip("fig07", "missing one of exp04/exp05/exp06 outputs")
        return False

    moea_pareto = np.asarray(
        json.loads(moea_json.read_text()).get("pareto_drought_chars") or []
    )
    lib_records = json.loads(lib_json.read_text())
    lib_feats = np.array([[r["mean_duration"], r["mean_avg_severity"]]
                          for r in lib_records], dtype=float)

    def _select(path: Path) -> np.ndarray:
        rec = json.loads(path.read_text())
        ids = set(rec.get("selected_trace_ids", []))
        return np.array(
            [[r["mean_duration"], r["mean_avg_severity"]]
             for r in lib_records if r.get("trace_id") in ids],
            dtype=float,
        )

    sets = [("MOEA-FIND", moea_pareto, COLORS["empirical"])]
    if lhs_json.exists():
        sets.append(("Library-LHS", _select(lhs_json), COLORS["lhs"]))
    if sobol_json.exists():
        sets.append(("Library-Sobol", _select(sobol_json), COLORS["sobol"]))

    fig = fig7_coverage_comparison(sets)
    _save(fig, FIGURES / "fig07_coverage_comparison.pdf")
    return True


def fig_si_hyperplane() -> bool:
    """SI-1 — hyperplane residual verification."""
    p2 = OUTPUTS / "exp01_analytic_2d" / "pareto.npz"
    p3 = OUTPUTS / "exp02_analytic_3d" / "pareto.npz"
    o2 = np.load(p2)["objs"] if p2.exists() else None
    o3 = np.load(p3)["objs"] if p3.exists() else None
    if o2 is None and o3 is None:
        _skip("figSI01", "neither 2D nor 3D Pareto available")
        return False
    fig = fig_si_hyperplane_check(o2, o3)
    _save(fig, FIGURES / "figSI01_hyperplane.pdf")
    return True


FACTORIES: Dict[str, Callable[[], bool]] = {
    "fig01": fig01_concept,
    "fig02_2d": fig02_2d_analytic,
    "fig02_3d": fig02_3d_analytic,
    "fig03": fig03_eps_nfe,
    "fig05": fig05_kirsch_pareto,
    "fig06": fig06_plausibility,
    "fig07": fig07_coverage,
    "figSI01": fig_si_hyperplane,
}


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--only", nargs="+", choices=sorted(FACTORIES.keys()),
                   help="Restrict to a subset of figures.")
    args = p.parse_args()

    keys = args.only or list(FACTORIES.keys())
    print(f"[10] assembling manuscript figures: {', '.join(keys)}")
    FIGURES.mkdir(parents=True, exist_ok=True)
    results: Dict[str, bool] = {}
    for k in keys:
        print(f"--- {k} ---")
        try:
            results[k] = FACTORIES[k]()
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR in {k}: {exc}")
            results[k] = False

    ok = sum(1 for v in results.values() if v)
    print(f"\n[10] {ok}/{len(results)} figures written to {FIGURES.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
