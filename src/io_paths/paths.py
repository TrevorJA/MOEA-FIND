"""Stage-aware output and figure path resolution.

Every workflow driver writes numerical artifacts under
``outputs/<stage>/<driver>/<slug>/`` and figures under
``figures/<stage>/<driver>/<slug>/``. The two trees mirror the
``workflows/<stage>/<driver>.py`` layout so a driver, its outputs, and
its figures share the same triple ``(stage, driver, slug)``.

Drivers should call :func:`stage_output_dir` and :func:`stage_figure_dir`
rather than hard-coding ``outputs/...`` or ``figures/...`` strings.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

# src/io_paths/paths.py -> parents[2] is the project root
# (parents[0]=io_paths, parents[1]=src, parents[2]=MOEA-FIND).
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS_ROOT = PROJECT_ROOT / "outputs"
FIGURES_ROOT = PROJECT_ROOT / "figures"


def stage_output_dir(
    stage: str,
    driver: str,
    slug: Optional[str] = None,
    *,
    create: bool = True,
) -> Path:
    """Return ``outputs/<stage>/<driver>[/<slug>]/``.

    Args:
        stage: Stage folder name, e.g. ``"01_analytic_validation"``.
        driver: Driver script stem, e.g. ``"analytic_2d"``.
        slug: Optional variant slug. Omit for stage-level shared dirs
            (e.g. calibration JSON written once per stage).
        create: ``True`` (default) to ``mkdir(parents=True, exist_ok=True)``.

    Returns:
        Absolute :class:`Path` under ``OUTPUTS_ROOT``.
    """
    p = OUTPUTS_ROOT / stage / driver
    if slug:
        p = p / slug
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p


def stage_figure_dir(
    stage: str,
    driver: str,
    slug: Optional[str] = None,
    *,
    create: bool = True,
) -> Path:
    """Return ``figures/<stage>/<driver>[/<slug>]/``.

    Companion to :func:`stage_output_dir` for plotting drivers.
    """
    p = FIGURES_ROOT / stage / driver
    if slug:
        p = p / slug
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p


def manuscript_figure_dir(kind: str = "main", *, create: bool = True) -> Path:
    """Return ``figures/main/`` or ``figures/supplementary/``.

    Only ``workflows/99_manuscript_figures/make_figures.py`` should
    write here; every other driver writes under
    :func:`stage_figure_dir`.
    """
    if kind not in {"main", "supplementary"}:
        raise ValueError(f"manuscript figure kind must be 'main' or 'supplementary', got {kind!r}")
    p = FIGURES_ROOT / kind
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p
