"""Self-describing variant slugs for MOEA-FIND output directories.

A slug encodes every parameter that distinguishes one experiment run from
another, in a filesystem-safe and human-scannable form. Slugs are the
final segment of an output path under ``outputs/expNN_<stage>/<slug>/``.

Format
------
``{stage}__{key1}={val1}__{key2}={val2}__...__s={seed}``

The ``__`` separator is unambiguous (an underscore is not a valid
separator since SSI timescales like ``ssi=3`` already contain digits).
The leading stage tag (``moea``, ``library``, ``subsample``, ``analytic``,
``pywrdrb``, ``discovery``) keeps mixed listings legible.

Numeric NFE values are rendered with ``k``/``M`` suffixes
(``500k``, ``2M``) for readability; the ``parse_slug`` helper round-trips
back to integers.

Stage helpers
-------------
- :func:`moea_slug`        — single-site or multi-site MOEA runs (stage 04, 05)
- :func:`library_slug`     — Kirsch library builds (stage 03)
- :func:`subsample_slug`   — LHS/Sobol subsampling (stage 03)
- :func:`analytic_slug`    — analytic validation (stage 01)
- :func:`pywrdrb_slug`     — Pywr-DRB re-evaluation (stage 06)
- :func:`build_slug`       — generic builder; takes ``(stage, **kwargs)``

All stages (including stage 04) build output-path slugs through these
helpers; the former ``make_variant_slug`` legacy format was retired in
the spring-cleaning reorg.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional


_KEY_VAL_SEP = "="
_FIELD_SEP = "__"


def _format_nfe(nfe: int) -> str:
    """Render an NFE count as ``500k`` / ``2M`` style for readability."""
    if nfe % 1_000_000 == 0:
        return f"{nfe // 1_000_000}M"
    if nfe >= 1_000_000:
        return f"{nfe / 1_000_000:.1f}M".rstrip("0").rstrip(".")
    if nfe % 1_000 == 0:
        return f"{nfe // 1_000}k"
    return str(nfe)


def _parse_nfe(token: str) -> int:
    """Inverse of :func:`_format_nfe`."""
    token = token.strip()
    if token.endswith("M"):
        return int(round(float(token[:-1]) * 1_000_000))
    if token.endswith("k"):
        return int(round(float(token[:-1]) * 1_000))
    return int(token)


def _format_value(value: Any) -> str:
    """Stringify a value for inclusion in a slug.

    Booleans become ``true``/``false``; ``None`` becomes ``none``;
    floats keep four significant digits; everything else is ``str()``.
    """
    if value is None:
        return "none"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _sanitize(token: str) -> str:
    """Replace filesystem-unsafe characters in a slug token."""
    return re.sub(r"[^A-Za-z0-9._\-+]", "-", token)


def build_slug(stage: str, **fields: Any) -> str:
    """Build a slug ``stage__k=v__k=v__...`` from arbitrary fields.

    The fields are emitted in the order they are passed as keyword
    arguments (Python preserves kwarg ordering since 3.7). Callers should
    therefore pass the most identifying fields first. Pass ``None`` to
    omit a field; ``False`` is included as ``false``.

    Args:
        stage: Stage tag (``"moea"``, ``"library"``, ``"subsample"``, ...).
        **fields: ``key=value`` parameters that describe the run.

    Returns:
        Slug string. Filesystem-safe.
    """
    parts = [_sanitize(stage)]
    for key, value in fields.items():
        if value is None:
            continue
        parts.append(
            f"{_sanitize(key)}{_KEY_VAL_SEP}{_sanitize(_format_value(value))}"
        )
    return _FIELD_SEP.join(parts)


def parse_slug(slug: str) -> Dict[str, Any]:
    """Parse a slug into a dict ``{'_stage': ..., 'key1': 'val1', ...}``.

    Values are returned as strings; callers are responsible for
    re-typing (e.g. ``int(parsed['T'])``). ``nfe`` tokens are
    auto-detected and converted via :func:`_parse_nfe`.

    Args:
        slug: A slug previously built by :func:`build_slug`.

    Returns:
        Dict with reserved key ``"_stage"`` plus one entry per field.
    """
    fields = slug.split(_FIELD_SEP)
    if not fields:
        return {}
    out: Dict[str, Any] = {"_stage": fields[0]}
    for f in fields[1:]:
        if _KEY_VAL_SEP not in f:
            continue
        k, v = f.split(_KEY_VAL_SEP, 1)
        if k == "nfe":
            try:
                out[k] = _parse_nfe(v)
                continue
            except ValueError:
                pass
        out[k] = v
    return out


# ---------------------------------------------------------------------------
# Stage-specific helpers — preferred entry points
# ---------------------------------------------------------------------------


def moea_slug(
    *,
    mode: str,
    n_years: int,
    nfe: int,
    seed: int,
    ssi: Optional[int] = None,
    metrics: Optional[str] = None,
    cons: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Slug for a single-site or multi-site MOEA run.

    Example: ``moea__mode=residual__T=20__nfe=500k__ssi=3__metrics=primary__cons=hydro__s=42``

    Args:
        mode: DV injection mode (``"index"`` / ``"residual"``).
        n_years: Synthetic trace length.
        nfe: Function evaluations.
        seed: Random seed.
        ssi: SSI accumulation timescale (1, 3, 6, 12). Optional.
        metrics: Metric-set preset name (``"primary"``, etc.). Optional.
        cons: Constraint mode (``"hydro"``, ``"dv-l2"``, ``"none"``).
            Optional.
        extra: Additional fields appended in sorted key order.

    Returns:
        Slug string.
    """
    fields: Dict[str, Any] = {
        "mode": mode,
        "T": n_years,
        "nfe": _format_nfe(nfe),
        "ssi": ssi,
        "metrics": metrics,
        "cons": cons,
    }
    if extra:
        for k in sorted(extra):
            fields[k] = extra[k]
    fields["s"] = seed
    return build_slug("moea", **fields)


def library_slug(
    *,
    n_traces: int,
    n_years: int,
    seed: int,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Slug for a Kirsch library build.

    Example: ``library__N=10000__T=20__s=42``

    Args:
        n_traces: Library size.
        n_years: Per-trace length in years.
        seed: Random seed.
        extra: Additional fields appended in sorted key order.
    """
    fields: Dict[str, Any] = {"N": n_traces, "T": n_years}
    if extra:
        for k in sorted(extra):
            fields[k] = extra[k]
    fields["s"] = seed
    return build_slug("library", **fields)


def subsample_slug(
    *,
    src: str,
    method: str,
    n_select: int,
    seed: int,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Slug for an LHS/Sobol subsample of a Kirsch library.

    Example: ``subsample__src=library_N10000_T20_s42__method=lhs__n=200__s=7``

    Args:
        src: Source library slug or short identifier.
        method: Subsampling method (``"lhs"``, ``"sobol"``, ``"random"``).
        n_select: Number of traces selected.
        seed: Random seed for the subsampler.
        extra: Additional fields appended in sorted key order.
    """
    fields: Dict[str, Any] = {"src": src, "method": method, "n": n_select}
    if extra:
        for k in sorted(extra):
            fields[k] = extra[k]
    fields["s"] = seed
    return build_slug("subsample", **fields)


def analytic_slug(
    *,
    k: int,
    nfe: int,
    seed: int,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Slug for an analytic validation run.

    Example: ``analytic__k=3__nfe=100k__s=0``

    Args:
        k: Number of objectives (dimensionality).
        nfe: Function evaluations.
        seed: Random seed.
        extra: Additional fields appended in sorted key order.
    """
    fields: Dict[str, Any] = {"k": k, "nfe": _format_nfe(nfe)}
    if extra:
        for key in sorted(extra):
            fields[key] = extra[key]
    fields["s"] = seed
    return build_slug("analytic", **fields)


def pywrdrb_slug(
    *,
    src: str,
    n_years: int,
    seed: int,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Slug for a Pywr-DRB re-evaluation run.

    Example: ``pywrdrb__src=moea_residual_T20_s42__T=20__s=42``

    Args:
        src: Identifier of the Pareto archive being re-evaluated (typically
            the upstream MOEA slug).
        n_years: Simulation length.
        seed: Random seed.
        extra: Additional fields appended in sorted key order.
    """
    fields: Dict[str, Any] = {"src": src, "T": n_years}
    if extra:
        for k in sorted(extra):
            fields[k] = extra[k]
    fields["s"] = seed
    return build_slug("pywrdrb", **fields)
