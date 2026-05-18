"""Stage 99 - Manuscript figure assembler.

Reads existing upstream stage outputs and regenerates the manuscript
figure set under ``figures/main/`` and ``figures/supplementary/``.
The script is idempotent and stateless: it regenerates every figure
whose upstream data exists and cleanly skips any figure whose upstream
data is missing.

Main-text figure inventory:

    Figure 1  fig01_param_vs_hazard_space.pdf   Main  §1.3     conceptual
    Figure 2  fig02_pipeline_schematic.pdf      Main  §2.6     Inkscape (manual)
    Figure 3  fig03_manhattan_construction.pdf  Main  §2.4     no upstream data
    Figure 4  fig04_dimension_sweep.pdf         Main  §3.1     01_analytic_validation/dimension_sweep
    Figure 5  fig05_cannonsville_hydrology.pdf  Main  §3.2     04_moea_find_single_site/run_moea_find/<slug>
    Figure 6  fig06_cannonsville_pareto.pdf     Main  §3.2     04_moea_find_single_site/run_moea_find/<slug>
                                                               + 03_kirsch_library/build_library
    Figure 7  fig07_gbt_hazard_discovery.pdf    Main  §3.3     06_pywrdrb_reeval/policy_reeval/<slug>

    figSI01_hyperplane.pdf                      SI-1           01_analytic_validation/analytic_2d, analytic_3d

Run:
    python workflows/99_manuscript_figures/make_figures.py
    python workflows/99_manuscript_figures/make_figures.py --only fig04
    python workflows/99_manuscript_figures/make_figures.py \
        --single-site-slug residual_T20_nfe200000_s42_constrained
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Callable, Dict

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

# Analytic figures run without SynHydro installed; figures that need real
# SSI calculations (fig05 rebuild, fig06) require SynHydro. Prefer the real
# module when it's importable so those paths work; only fall back to a
# MagicMock stub when SynHydro is genuinely missing.
try:  # noqa: E402
    import synhydro  # noqa: F401
    import synhydro.droughts  # noqa: F401
    import synhydro.droughts.ssi  # noqa: F401
except Exception:
    from unittest.mock import MagicMock  # noqa: E402
    _stub = MagicMock()
    sys.modules.setdefault("synhydro", _stub)
    sys.modules.setdefault("synhydro.droughts", _stub.droughts)
    sys.modules.setdefault("synhydro.droughts.ssi", _stub.droughts.ssi)

import matplotlib.pyplot as plt  # noqa: E402

from src.io_paths.paths import (  # noqa: E402
    OUTPUTS_ROOT,
    manuscript_figure_dir,
    stage_output_dir,
)
from src.plotting.analytic import (  # noqa: E402
    fig3_manhattan_construction,
    fig4_dimension_sweep,
    fig_si_hyperplane_check,
)
from src.plotting.architecture import fig3_wrapper_schematic  # noqa: E402

OUTPUTS = OUTPUTS_ROOT
FIGURES_MAIN = manuscript_figure_dir("main")
FIGURES_SI = manuscript_figure_dir("supplementary")

# Stage roots used by the figure factories below. ``stage_output_dir``
# without a slug returns the driver-level dir (a stage-shared parent).
DIMENSION_SWEEP_DIR: Path = stage_output_dir(
    "01_analytic_validation", "dimension_sweep", create=False,
)
ANALYTIC_2D_DIR: Path = stage_output_dir(
    "01_analytic_validation", "analytic_2d", create=False,
)
ANALYTIC_3D_DIR: Path = stage_output_dir(
    "01_analytic_validation", "analytic_3d", create=False,
)
KIRSCH_LIBRARY_DIR: Path = stage_output_dir(
    "03_kirsch_library", "build_library", create=False,
)
# Production library slug. Override with --library-slug if a different
# build was used. fig06 needs this to find characteristics.npz.
KIRSCH_LIBRARY_SLUG: str = "n10000_t20_ssi3_s42"
# Analytic compute slugs (figSI01 reads pareto.npz from these).
ANALYTIC_2D_SLUG: str = "k2_nfe100000_eps0.060_s42"
ANALYTIC_3D_SLUG: str = "k3_nfe50000_eps0.150_s42"

# Single-site MOEA-FIND variant dir + matching pywrdrb re-eval dir.
# Resolved at main() time from ``--single-site-slug``.
SINGLE_SITE_DIR: Path = stage_output_dir(
    "04_moea_find_single_site", "run_moea_find", create=False,
)
POLICY_REEVAL_DIR: Path = stage_output_dir(
    "06_pywrdrb_reeval", "policy_reeval", create=False,
)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(PROJECT_ROOT)}")


def _skip(name: str, reason: str) -> None:
    print(f"  SKIP {name}: {reason}")


# =============================================================================
# Generatable now: figures 1, 3, 4 and SI-1
# =============================================================================

def fig01_param_vs_hazard() -> bool:
    """Figure 1 (§1.3) - parameter space versus drought hazard space.

    Deferred to the final polish pass: the conceptual narrative for this
    figure (qualitative: parameter space space-fills, hazard space does
    not) is hard to communicate faithfully without real Kirsch-Nowak
    library output, and the previous synthetic fallback misrepresented
    the generator. A blank placeholder is emitted so downstream LaTeX
    compiles; the real panel will land after Figures 2-4 stabilise.
    """
    fig = plt.figure(figsize=(7.0, 3.0))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(
        0.5, 0.6,
        "Figure 1 placeholder",
        ha="center", va="center", fontsize=14, weight="bold",
        color="#555555",
    )
    ax.text(
        0.5, 0.38,
        "parameter space vs drought hazard space - deferred\n"
        "(will use real Kirsch-Nowak library projection)",
        ha="center", va="center", fontsize=9, style="italic",
        color="#777777",
    )
    fig.patch.set_edgecolor("#cccccc")
    fig.patch.set_linewidth(1.0)
    _save(fig, FIGURES_MAIN / "fig01_param_vs_hazard_space.pdf")
    return True


def fig03_manhattan() -> bool:
    """Figure 3 (§2.4) - Manhattan-distance auxiliary objective construction."""
    fig = fig3_manhattan_construction()
    _save(fig, FIGURES_MAIN / "fig03_manhattan_construction.pdf")
    return True


def fig04_dim_sweep() -> bool:
    """Figure 4 (§3.1) - interior-filling coverage at K = 2 through K = 6."""
    if not DIMENSION_SWEEP_DIR.exists():
        _skip("fig04", f"missing {DIMENSION_SWEEP_DIR.relative_to(PROJECT_ROOT)}")
        return False
    fig = fig4_dimension_sweep(DIMENSION_SWEEP_DIR)
    _save(fig, FIGURES_MAIN / "fig04_dimension_sweep.pdf")
    return True


# =============================================================================
# Figure 2 pipeline schematic (programmatic placeholder; Inkscape later)
# =============================================================================
def fig02_pipeline() -> bool:
    """Figure 2 (§2.6) - MOEA-FIND algorithmic pipeline.

    Programmatic matplotlib schematic rendered from
    :func:`src.plotting.architecture.fig3_wrapper_schematic` as a
    self-generating placeholder. Will be replaced by an Inkscape-authored
    SVG in the final polish pass.
    """
    fig = fig3_wrapper_schematic()
    _save(fig, FIGURES_MAIN / "fig02_pipeline_schematic.pdf")
    return True


# =============================================================================
# Manuscript figures 5, 6, 7 - sourced from stage 04 / stage 06 outputs
# =============================================================================
#
# fig05 and fig06 are thin re-packagers: they copy per-variant PDFs
# produced by the stage 04 ``run_moea_find`` driver into
# ``figures/main/`` with manuscript-ordered names. No new plotting is
# done here; if the upstream PDFs are missing the factory returns False
# with a helpful SKIP reason. fig07 packages the satisficing-map PDF
# from stage 06 ``policy_reeval``.

def _copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    print(f"  wrote {dst.relative_to(PROJECT_ROOT)} "
          f"(from {src.relative_to(PROJECT_ROOT)})")


def fig05_cannonsville_hydrology() -> bool:
    """Figure 5 (§3.2) - Cannonsville hydrology statistical-agreement panel.

    Rebuilds the 2x2 combined panel (ACF, FDC, seasonal mean, seasonal
    std) from the synthetic Pareto traces stored in
    ``SINGLE_SITE_DIR/results.json`` and freshly resampled historical
    blocks at the same block length. If a prebuilt ``fig05_hydrology.pdf``
    already sits in ``SINGLE_SITE_DIR/figures``, that is copied instead.
    """
    import json
    src_dir = SINGLE_SITE_DIR / "figures"
    results_path = SINGLE_SITE_DIR / "results.json"

    # Fast path: the upstream driver already produced the combined PDF.
    prebuilt = src_dir / "fig05_hydrology.pdf"
    if prebuilt.exists():
        _copy(prebuilt, FIGURES_MAIN / "fig05_hydrology.pdf")
        return True

    if not results_path.exists():
        _skip("fig05", f"missing {results_path.relative_to(PROJECT_ROOT)} "
                        "(run stage 04 run_moea_find for the configured "
                        "variant first)")
        return False

    # Rebuild the combined panel from on-disk traces + fresh historical blocks.
    try:
        from src.plotting.trace_diagnostics import plot_hydrology_panels
        from src.hydrology.historical_blocks import (
            resample_historical_blocks,
            resample_historical_blocks_2d,
        )
        from src.experiment import prepare_data
    except Exception as exc:  # noqa: BLE001
        _skip("fig05", f"imports for rebuild failed: {exc}")
        return False

    results = json.loads(results_path.read_text())
    syn_1d_raw = results.get("pareto_traces_1d") or []
    syn_2d_raw = results.get("pareto_traces_2d") or []
    if not syn_1d_raw or not syn_2d_raw:
        _skip("fig05", "results.json has no pareto_traces_1d/_2d")
        return False
    n_years = int(results.get("n_years_out", len(syn_2d_raw[0])))
    syn_1d = [np.asarray(t) for t in syn_1d_raw]
    syn_2d = [np.asarray(t) for t in syn_2d_raw]

    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    monthly_2d, monthly_1d = prepare_data(cache_dir)

    hist_1d = resample_historical_blocks(
        monthly_1d, T_years=n_years, stride=1,
    )
    hist_2d = resample_historical_blocks_2d(
        monthly_2d, T_years=n_years, stride=1,
    )
    fig, _ = plot_hydrology_panels(syn_1d, hist_1d, syn_2d, hist_2d)
    _save(fig, FIGURES_MAIN / "fig05_hydrology.pdf")
    return True


def fig06_cannonsville_pareto() -> bool:
    """Figure 6 (§3.2) - Cannonsville drought-space composite.

    Composite figure: 2D scatter + marginals (three populations), 3D
    scatter (three populations), and four sample monthly-flow traces
    (two historical blocks at the drought extremes + the nearest
    MOEA-FIND Pareto traces in 2D drought space), with SSI-based
    drought events shaded.

    Fast path: copies ``SINGLE_SITE_DIR/figures/fig06_composite.pdf``
    if present (emitted by run_moea_find). Otherwise rebuilds from
    ``SINGLE_SITE_DIR/results.json``,
    ``KIRSCH_LIBRARY_DIR/characteristics.npz``, and the cached USGS
    historical series.
    """
    import json
    src_dir = SINGLE_SITE_DIR / "figures"
    results_path = SINGLE_SITE_DIR / "results.json"
    prebuilt = src_dir / "fig06_composite.pdf"
    if prebuilt.exists():
        _copy(prebuilt, FIGURES_MAIN / "fig06_cannonsville_pareto.pdf")
        return True

    if not results_path.exists():
        _skip("fig06", f"missing {results_path.relative_to(PROJECT_ROOT)} "
                        "(run stage 04 run_moea_find for the configured "
                        "variant first)")
        return False

    kirsch_npz = KIRSCH_LIBRARY_DIR / KIRSCH_LIBRARY_SLUG / "characteristics.npz"
    if not kirsch_npz.exists():
        _skip("fig06", f"missing {kirsch_npz.relative_to(PROJECT_ROOT)} "
                        "(run stage 03 build_library first)")
        return False

    try:
        from src.plotting.drought_space import plot_fig6_composite
        from src.hydrology.historical_blocks import compute_historical_block_chars
        from src.experiment import prepare_data
        from src.metrics.objectives import make_ssi_calculator, flows_to_series
    except Exception as exc:  # noqa: BLE001
        _skip("fig06", f"imports for rebuild failed: {exc}")
        return False

    results = json.loads(results_path.read_text())
    pareto_metrics = np.asarray(results["drought_metrics"], dtype=float)
    pareto_traces_1d_raw = results.get("pareto_traces_1d") or []
    if not pareto_traces_1d_raw:
        _skip("fig06", "results.json has no pareto_traces_1d")
        return False
    pareto_traces_1d = [np.asarray(t) for t in pareto_traces_1d_raw]
    anti_ideal = np.asarray(results.get("anti_ideal", []), dtype=float)
    objective_keys = list(results.get("objective_keys") or [])
    n_years = int(results.get("n_years_out", 20))
    ssi_timescale = int(results.get("ssi_timescale", 3))

    # Kirsch library ---------------------------------------------------------
    kdata = dict(np.load(kirsch_npz))
    kirsch_all_keys = [str(k) for k in kdata["all_keys"]]
    kirsch_all_values = np.asarray(kdata["all_values"], dtype=float)

    def _kirsch_col(name):
        if name not in kirsch_all_keys:
            return None
        return kirsch_all_values[:, kirsch_all_keys.index(name)]

    # Project to objective-key order used by Borg (typically 3 keys).
    kirsch_by_objective = []
    for k in objective_keys:
        col = _kirsch_col(k)
        if col is None:
            print(f"  note: Kirsch library lacks '{k}'; dropping axis from 3D")
            kirsch_by_objective = []
            break
        kirsch_by_objective.append(col)
    if kirsch_by_objective:
        kirsch_full = np.column_stack(kirsch_by_objective)
    else:
        kirsch_full = np.empty((0, len(objective_keys)))
    # 2D kirsch projection is the first two objective-key columns.
    if kirsch_full.shape[1] >= 2:
        kirsch_2d = kirsch_full[:, :2]
    else:
        kirsch_2d = np.empty((0, 2))
    kirsch_3d = kirsch_full if kirsch_full.shape[1] >= 3 else None

    # Historical record + SSI calculator ------------------------------------
    cache_dir = PROJECT_ROOT / "outputs" / "data_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    _, monthly_1d = prepare_data(cache_dir)

    ssi_calc = make_ssi_calculator(timescale=ssi_timescale)
    ssi_calc.fit(flows_to_series(monthly_1d))

    # Historical block drought characteristics, same objective-key order.
    hist_chars = compute_historical_block_chars(
        monthly_1d, T_years=n_years, ssi_calc=ssi_calc,
        objective_keys=objective_keys, stride=1,
    )

    # 1D historical block traces at stride=1 (same ordering as hist_chars).
    from src.hydrology.historical_blocks import resample_historical_blocks
    hist_blocks_1d = resample_historical_blocks(
        monthly_1d, T_years=n_years, stride=1,
    )

    # Recover the water-year calendar so panels (c)/(d) label real dates.
    # prepare_data aligns monthly_1d to the first October of the daily
    # USGS cache; the daily cache stores a datetime index from which we
    # can read the first year cheaply.
    import pandas as pd
    try:
        daily_csv = next(cache_dir.glob("usgs_*_daily.csv"))
        daily_df = pd.read_csv(daily_csv, index_col=0, parse_dates=True,
                                nrows=500)
        first_oct = daily_df.index[daily_df.index.month == 10][0]
        hist_wy_start = int(first_oct.year)
    except Exception:
        hist_wy_start = None
    if hist_wy_start is not None:
        # Block i covers water years [start_year+i, start_year+i+T_years).
        start_years = [hist_wy_start + i for i in range(len(hist_blocks_1d))]
    else:
        start_years = None

    fig = plot_fig6_composite(
        pareto_metrics=pareto_metrics,
        pareto_traces_1d=pareto_traces_1d,
        kirsch_2d=kirsch_2d,
        kirsch_3d=kirsch_3d,
        hist_block_chars_2d=hist_chars[:, :2],
        hist_block_chars_3d=hist_chars if hist_chars.shape[1] >= 3 else None,
        hist_blocks_1d=hist_blocks_1d,
        ssi_calc=ssi_calc,
        anti_ideal=anti_ideal if anti_ideal.size else None,
        historical_point_2d=None,
        historical_point_3d=None,
        objective_labels_2d=(
            objective_keys[0] if objective_keys else "D_1",
            objective_keys[1] if len(objective_keys) > 1 else "D_2",
        ),
        objective_labels_3d=tuple(
            objective_keys[:3] + [""] * max(0, 3 - len(objective_keys))
        )[:3],
        historical_start_years=start_years,
    )
    _save(fig, FIGURES_MAIN / "fig06_cannonsville_pareto.pdf")
    return True


def fig07_gbt_hazard_discovery() -> bool:
    """Figure 7 (§3.3) - DRB scenario discovery (satisficing map).

    Repackages the satisficing classification figure emitted by
    ``src/plotting/07_scenario_discovery/scenario_discovery_plots.py``
    (under ``figures/07_scenario_discovery/scenario_discovery_plots/<slug>/``)
    as the manuscript-ordered ``fig07_gbt_hazard_discovery.pdf``.
    """
    from src.io_paths.paths import stage_figure_dir

    sd_dir = stage_figure_dir(
        "07_scenario_discovery", "scenario_discovery_plots",
        SINGLE_SITE_DIR.name,  # slug matches stage 04 single-site slug
        create=False,
    )
    candidates = [
        sd_dir / "fig09_satisficing_map.pdf",
        sd_dir / "fig09_satisficing_map_gbt.pdf",
    ]
    for prebuilt in candidates:
        if prebuilt.exists():
            _copy(prebuilt, FIGURES_MAIN / "fig07_gbt_hazard_discovery.pdf")
            return True
    _skip("fig07", f"none of {[p.relative_to(PROJECT_ROOT) for p in candidates]} "
                    "exist (run stage 07 scenario_discovery_plots first)")
    return False


# =============================================================================
# Supporting Information figures
# =============================================================================
def fig_si01_hyperplane() -> bool:
    """Figure SI-1 - hyperplane residual histogram."""
    p2 = ANALYTIC_2D_DIR / ANALYTIC_2D_SLUG / "pareto.npz"
    p3 = ANALYTIC_3D_DIR / ANALYTIC_3D_SLUG / "pareto.npz"
    o2 = np.load(p2)["objs"] if p2.exists() else None
    o3 = np.load(p3)["objs"] if p3.exists() else None
    if o2 is None and o3 is None:
        _skip("figSI01", "neither 2D nor 3D Pareto available")
        return False
    fig = fig_si_hyperplane_check(o2, o3)
    _save(fig, FIGURES_SI / "figSI01_hyperplane.pdf")
    return True


FACTORIES: Dict[str, Callable[[], bool]] = {
    "fig01": fig01_param_vs_hazard,
    "fig02": fig02_pipeline,
    "fig03": fig03_manhattan,
    "fig04": fig04_dim_sweep,
    "fig05": fig05_cannonsville_hydrology,
    "fig06": fig06_cannonsville_pareto,
    "fig07": fig07_gbt_hazard_discovery,
    "figSI01": fig_si01_hyperplane,
}


def main():
    global SINGLE_SITE_DIR, POLICY_REEVAL_DIR
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--only", nargs="+", choices=sorted(FACTORIES.keys()),
                   help="Restrict to a subset of figures.")
    p.add_argument(
        "--single-site-slug", type=str, default=None,
        help="Variant slug under "
             "outputs/04_moea_find_single_site/run_moea_find/<slug>/. "
             "The same slug is used to locate the matching stage 06 "
             "outputs/06_pywrdrb_reeval/policy_reeval/<slug>/. "
             "Required for fig05, fig06, and fig07.",
    )
    args = p.parse_args()

    if args.single_site_slug:
        SINGLE_SITE_DIR = stage_output_dir(
            "04_moea_find_single_site", "run_moea_find",
            args.single_site_slug, create=False,
        )
        POLICY_REEVAL_DIR = stage_output_dir(
            "06_pywrdrb_reeval", "policy_reeval",
            args.single_site_slug, create=False,
        )
        print(f"[99] single-site slug: {args.single_site_slug}")
        print(f"[99] single-site dir:  "
              f"{SINGLE_SITE_DIR.relative_to(PROJECT_ROOT)}")
        print(f"[99] policy-reeval:    "
              f"{POLICY_REEVAL_DIR.relative_to(PROJECT_ROOT)}")
    else:
        print("[99] no --single-site-slug provided; "
              "fig05/fig06/fig07 will SKIP unless their default driver "
              "directory contains results.")

    keys = args.only or list(FACTORIES.keys())
    print(f"[99] assembling manuscript figures: {', '.join(keys)}")
    print(f"[99] main:          {FIGURES_MAIN.relative_to(PROJECT_ROOT)}")
    print(f"[99] supplementary: {FIGURES_SI.relative_to(PROJECT_ROOT)}")
    results: Dict[str, bool] = {}
    for k in keys:
        print(f"--- {k} ---")
        try:
            results[k] = FACTORIES[k]()
        except NotImplementedError as exc:
            print(f"  PENDING {k}: {exc}")
            results[k] = False
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR in {k}: {exc}")
            import traceback
            traceback.print_exc()
            results[k] = False

    ok = sum(1 for v in results.values() if v)
    print(f"\n[99] {ok}/{len(results)} figures written under "
          f"{FIGURES_MAIN.parent.relative_to(PROJECT_ROOT)}/")


if __name__ == "__main__":
    main()
