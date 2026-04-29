"""Stage 08 — Global sensitivity analysis of NYC reservoir outcomes.

Consumes the Stage-06 metric bank (Y) and the upstream MOEA-FIND
realization archive's drought characteristics (X), then runs three
sample-free SA methods — Delta moment-independent (default), PAWN,
RBD-FAST — on every chosen NYC outcome. For each (outcome, method),
emits headline indices, bootstrap CIs, rank-stability, and the inputs
needed for cross-method and cross-outcome diagnostics.

Factor space = the optimized MOEA-FIND objective axes (read from
``results.json::objective_keys`` by default; overrideable with
``--metric-set``). Non-optimized drought characteristics are excluded
because they have no space-filling guarantee across the realization
sample.

Compute-only driver: writes numerical artifacts under
``outputs/08_nyc_sensitivity/run_sa/<slug>/``. Figures are produced by
the paired plotting driver
``workflows/08_nyc_sensitivity/plots/run_sa.py``.

Outputs (``outputs/08_nyc_sensitivity/run_sa/<slug>/``):

    config.json
    results/
        indices_<method>.parquet            (long-form: outcome × factor × cols)
        bootstrap_<method>.parquet
        rank_stability_<method>.parquet
        rank_stability_summary_<method>.parquet
        convergence_<method>.parquet
        cross_method_rank_corr.parquet      (per outcome, method × method)
        cross_outcome_rank_corr.parquet     (per method, outcome × outcome)
        selection_log.json

Usage:
    python workflows/08_nyc_sensitivity/run_sa.py \\
        --bank outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/metric_bank.parquet \\
        --chars outputs/04_moea_find_single_site/run_moea_find/<src_slug>/results.json \\
        --config workflows/08_nyc_sensitivity/configs/all_methods.yaml
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from src.drought_metrics import PRESETS, REGISTRY, resolve_metric_set, metric_names  # noqa: E402
from src.io import load_metric_bank, load_pareto_chars, save_experiment_config  # noqa: E402
from src.paths import stage_output_dir  # noqa: E402
from src.sensitivity import (  # noqa: E402
    HEADLINE_INDEX,
    METHODS,
    apply_method_selection_criterion,
    bootstrap_indices,
    bootstrap_rank_stability,
    convergence_curve,
    cross_method_rank_corr,
    cross_outcome_rank_corr,
    drop_low_variance_factors,
    resolve_method,
)
from src.slugs import build_slug  # noqa: E402

STAGE = "08_nyc_sensitivity"
DRIVER = "run_sa"

DEFAULT_METHODS = ("delta", "pawn", "rbd_fast")
DEFAULT_OUTCOMES = (
    "nyc_min_storage_frac",
    "nyc_drawdown_days_below_0.25",
    "montague_flow_reliability",
    "montague_flow_vulnerability",
)
DEFAULT_LOG_TRANSFORM = (
    "nyc_drawdown_days_below_0.25",
    "montague_flow_vulnerability",
)
DEFAULT_CONVERGENCE_SIZES = (200, 500, 1000)


_YAML_TO_CLI = {
    "methods": "methods",
    "outcomes": "outcomes",
    "metric_set": "metric_set",
    "n_bootstrap": "n_bootstrap",
    "n_convergence_replicates": "n_convergence_replicates",
    "convergence_sizes": "convergence_sizes",
    "seed": "seed",
    "cv_drop_threshold": "cv_drop_threshold",
    "log_transform_outcomes": "log_transform_outcomes",
    "rank_spearman_threshold": "rank_spearman_threshold",
    "cross_method_threshold": "cross_method_threshold",
    "delta_num_resamples": "delta_num_resamples",
    "pawn_S": "pawn_S",
    "rbd_fast_M": "rbd_fast_M",
    "rbd_fast_num_resamples": "rbd_fast_num_resamples",
}


def _load_yaml_config(config_path: Path) -> dict:
    """Load a YAML preset into argparse-compatible defaults.

    Unknown keys are warned but not fatal so a typo doesn't silently
    misconfigure a run.
    """
    import yaml

    raw = yaml.safe_load(Path(config_path).read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{config_path}: top-level YAML must be a mapping")
    overrides: dict = {}
    for k, v in raw.items():
        if k in _YAML_TO_CLI:
            overrides[_YAML_TO_CLI[k]] = v
        else:
            print(f"[08] WARN: unknown YAML key {k!r} in {config_path}; ignored")
    return overrides


# ---------------------------------------------------------------------------
# Factor + outcome resolution
# ---------------------------------------------------------------------------


def _resolve_factor_set(
    chars_payload: dict,
    metric_set_arg: Optional[str],
) -> Tuple[List[str], str]:
    """Resolve the SA factor set + slug-friendly metric-set tag.

    The factor set is the MOEA-FIND objective axes. By default these are
    read from ``chars_payload['objective_keys']`` (the upstream archive's
    record of what was optimized). The user can override with
    ``--metric-set <preset|comma-list>``; the override is validated
    against ``src.drought_metrics.REGISTRY``.

    Returns:
        Tuple of ``(factor_names, metric_set_tag)``. ``metric_set_tag``
        is a slug-safe identifier of the factor set: a preset name when
        the input is a registered preset, else a short SHA-1 prefix of
        the comma-joined factor names.
    """
    if metric_set_arg is None:
        obj_keys = list(chars_payload.get("objective_keys") or [])
        if not obj_keys:
            raise SystemExit(
                "[08] no --metric-set supplied and chars JSON has no "
                "'objective_keys'. Supply --metric-set explicitly or "
                "regenerate the upstream archive with objective_keys."
            )
        # Drop legacy Manhattan-norm slot if it sneaks in. The Manhattan
        # objective is not a registered metric in REGISTRY, so any name
        # outside REGISTRY is a candidate for dropping.
        dropped = [k for k in obj_keys if k not in REGISTRY]
        if dropped:
            print(f"[08] dropping non-registry objective_keys entries: {dropped}")
        factor_names = [k for k in obj_keys if k in REGISTRY]
        if not factor_names:
            raise SystemExit(
                f"[08] objective_keys {obj_keys} contains no entries in "
                "drought_metrics.REGISTRY. Cannot infer factor set; pass "
                "--metric-set explicitly."
            )
        # Tag from the upstream metric_set field if present, else hash.
        ms_field = chars_payload.get("metric_set")
        if isinstance(ms_field, str) and ms_field in PRESETS:
            tag = ms_field
        else:
            tag = "h" + hashlib.sha1(
                ",".join(factor_names).encode("utf-8")
            ).hexdigest()[:6]
        return factor_names, tag

    # Explicit override.
    if metric_set_arg in PRESETS:
        metric_set = resolve_metric_set(metric_set_arg)
        return list(metric_names(metric_set)), metric_set_arg
    # Comma-list of metric names
    names = [n.strip() for n in metric_set_arg.split(",") if n.strip()]
    metric_set = resolve_metric_set(names)
    factor_names = list(metric_names(metric_set))
    tag = "h" + hashlib.sha1(
        ",".join(factor_names).encode("utf-8")
    ).hexdigest()[:6]
    return factor_names, tag


def _align_xy(
    chars_df: pd.DataFrame,
    bank_df: pd.DataFrame,
    factor_names: Sequence[str],
    outcome: str,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Inner-join chars and bank by realization id; build (X, Y, kept_ids)."""
    common = chars_df.index.intersection(bank_df.index)
    if not len(common):
        raise SystemExit(
            f"[08] no realization ids in common between chars "
            f"({len(chars_df)}) and metric bank ({len(bank_df)})."
        )
    X = chars_df.loc[common, list(factor_names)].astype(float)
    Y = pd.to_numeric(bank_df.loc[common, outcome], errors="coerce")
    finite = X.notna().all(axis=1) & Y.notna()
    return (
        X.loc[finite].values,
        Y.loc[finite].values.astype(float),
        list(common[finite]),
    )


# ---------------------------------------------------------------------------
# Method invocation
# ---------------------------------------------------------------------------


def _method_kwargs(method_name: str, args) -> Dict:
    """Per-method extra kwargs harvested from CLI / YAML."""
    if method_name == "delta":
        return {
            "num_resamples": args.delta_num_resamples,
            "conf_level": 0.95,
        }
    if method_name == "pawn":
        return {
            "S": args.pawn_S,
            "n_bootstrap": 0,  # bootstrap CI handled by diagnostic layer
        }
    if method_name == "rbd_fast":
        return {
            "M": args.rbd_fast_M,
            "num_resamples": args.rbd_fast_num_resamples,
            "conf_level": 0.95,
        }
    return {}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None,
                     help="YAML preset under workflows/08_nyc_sensitivity/configs/.")
    pre_args, _ = pre.parse_known_args()
    yaml_defaults = _load_yaml_config(pre_args.config) if pre_args.config else {}

    p = argparse.ArgumentParser(
        parents=[pre],
        description="Stage 08 — global sensitivity analysis of NYC outcomes.",
    )
    p.add_argument("--bank", type=Path, required=True,
                   help="Stage-06 metric_bank.parquet under "
                        "outputs/06_pywrdrb_reeval/policy_reeval/<src_slug>/results/.")
    p.add_argument("--chars", type=Path, required=True,
                   help="Upstream MOEA-FIND results.json under "
                        "outputs/04_moea_find_single_site/run_moea_find/<src_slug>/. "
                        "Must contain 'pareto_chars' and 'objective_keys'.")
    p.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS),
                   choices=sorted(METHODS.keys()),
                   help="SA methods to run. Default: all three.")
    p.add_argument("--outcomes", nargs="+", default=list(DEFAULT_OUTCOMES),
                   help="Metric-bank columns to analyze.")
    p.add_argument("--metric-set", default=None,
                   help="Preset name (primary, extreme_event, trace_fdc, "
                        "legacy) OR comma-separated metric names. If omitted, "
                        "the factor set is read from the upstream chars JSON's "
                        "objective_keys.")
    p.add_argument("--n-bootstrap", type=int, default=1000,
                   help="Bootstrap replicates for CI / rank-stability.")
    p.add_argument("--n-convergence-replicates", type=int, default=200,
                   help="Random sub-sample replicates per convergence size.")
    p.add_argument("--convergence-sizes", nargs="+", type=int,
                   default=list(DEFAULT_CONVERGENCE_SIZES))
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--cv-drop-threshold", type=float, default=0.05,
                   help="Drop factors whose |std/mean| is below this CV.")
    p.add_argument("--log-transform-outcomes", nargs="*",
                   default=list(DEFAULT_LOG_TRANSFORM),
                   help="Outcomes to additionally analyze on log1p scale.")
    p.add_argument("--rank-spearman-threshold", type=float, default=0.8)
    p.add_argument("--cross-method-threshold", type=float, default=0.7)
    p.add_argument("--delta-num-resamples", type=int, default=100)
    p.add_argument("--pawn-S", type=int, default=10)
    p.add_argument("--rbd-fast-M", type=int, default=10)
    p.add_argument("--rbd-fast-num-resamples", type=int, default=100)
    if yaml_defaults:
        p.set_defaults(**yaml_defaults)
    args = p.parse_args()

    # --- Load X (chars) and Y (metric bank) ---
    print(f"[08] loading chars from {args.chars}")
    chars_payload = json.loads(Path(args.chars).read_text())
    chars_df = load_pareto_chars(args.chars)

    print(f"[08] loading metric bank from {args.bank}")
    bank_df = load_metric_bank(args.bank)

    # --- Resolve factor set (the optimized MOEA-FIND axes) ---
    factor_names, metric_set_tag = _resolve_factor_set(
        chars_payload, args.metric_set
    )
    print(f"[08] factor set ({metric_set_tag}): {factor_names}")

    # Verify every requested factor is in the chars frame.
    missing = [f for f in factor_names if f not in chars_df.columns]
    if missing:
        raise SystemExit(
            f"[08] factor names {missing} are not columns of pareto_chars. "
            f"Available: {list(chars_df.columns)}"
        )

    # --- Drop low-variance factors (logged + persisted) ---
    factor_view = chars_df[list(factor_names)]
    factor_view_kept, dropped_factors = drop_low_variance_factors(
        factor_view, cv_threshold=args.cv_drop_threshold,
    )
    if dropped_factors:
        print(f"[08] dropping low-variance factors (CV<"
              f"{args.cv_drop_threshold}): {dropped_factors}")
    factor_names = list(factor_view_kept.columns)
    if not factor_names:
        raise SystemExit("[08] all factors dropped; nothing to analyze.")

    # --- Validate outcome list against the bank ---
    missing_outcomes = [o for o in args.outcomes if o not in bank_df.columns]
    if missing_outcomes:
        raise SystemExit(
            f"[08] outcomes {missing_outcomes} not in metric bank. "
            f"Available: {sorted(bank_df.columns.tolist())}"
        )

    # --- Slug + output paths ---
    upstream_slug = chars_payload.get("variant") or Path(args.chars).parent.name
    slug = build_slug(
        "sa",
        src=upstream_slug,
        metric_set=metric_set_tag,
        methods="-".join(sorted(args.methods)),
        n_outcomes=len(args.outcomes),
        s=args.seed,
    )
    out = stage_output_dir(STAGE, DRIVER, slug)
    results_dir = out / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    print(f"[08] variant: {slug}")
    print(f"[08] output: {out}")

    # --- Config dump ---
    save_experiment_config(out, {
        "script": "workflows/08_nyc_sensitivity/run_sa.py",
        "stage": STAGE,
        "driver": DRIVER,
        "variant": slug,
        "bank": str(args.bank),
        "chars": str(args.chars),
        "methods": list(args.methods),
        "outcomes": list(args.outcomes),
        "metric_set_tag": metric_set_tag,
        "factor_names": list(factor_names),
        "dropped_factors": dropped_factors,
        "n_bootstrap": args.n_bootstrap,
        "n_convergence_replicates": args.n_convergence_replicates,
        "convergence_sizes": list(args.convergence_sizes),
        "seed": args.seed,
        "cv_drop_threshold": args.cv_drop_threshold,
        "log_transform_outcomes": list(args.log_transform_outcomes),
        "selection": {
            "rank_spearman_threshold": args.rank_spearman_threshold,
            "cross_method_threshold": args.cross_method_threshold,
        },
    })

    # --- Build per-(outcome, method) result tables ---
    # Outcome list, including log1p variants.
    full_outcome_list: List[Tuple[str, str, str]] = []  # (label, raw_col, transform)
    for oc in args.outcomes:
        full_outcome_list.append((oc, oc, "raw"))
        if oc in args.log_transform_outcomes:
            full_outcome_list.append((f"{oc}_log1p", oc, "log1p"))

    # method_results[method][outcome_label] = per-factor DataFrame
    method_results: Dict[str, Dict[str, pd.DataFrame]] = {
        m: {} for m in args.methods
    }
    bootstrap_results: Dict[str, Dict[str, pd.DataFrame]] = {
        m: {} for m in args.methods
    }
    rank_stab_results: Dict[str, Dict[str, pd.DataFrame]] = {
        m: {} for m in args.methods
    }
    convergence_results: Dict[str, Dict[str, pd.DataFrame]] = {
        m: {} for m in args.methods
    }

    realization_counts: Dict[str, int] = {}

    t0 = time.time()
    for outcome_label, raw_col, transform in full_outcome_list:
        X_arr, Y_arr, kept_ids = _align_xy(
            chars_df, bank_df, factor_names, raw_col,
        )
        if transform == "log1p":
            Y_arr = np.log1p(Y_arr)
        realization_counts[outcome_label] = int(X_arr.shape[0])
        if X_arr.shape[0] < 50:
            print(f"[08] outcome {outcome_label}: only {X_arr.shape[0]} "
                  f"realizations after NaN drop; SA will be unreliable but "
                  f"still computed.")

        for method_name in args.methods:
            method = resolve_method(method_name)
            kw = _method_kwargs(method_name, args)
            print(f"[08] {outcome_label} / {method_name}: computing indices "
                  f"on n={X_arr.shape[0]}...")

            indices_df = method.compute(X_arr, Y_arr, factor_names, seed=args.seed, **kw)
            method_results[method_name][outcome_label] = indices_df

            boot_df = bootstrap_indices(
                method, X_arr, Y_arr, factor_names,
                n_bootstrap=args.n_bootstrap, seed=args.seed,
                method_kwargs=kw,
            )
            bootstrap_results[method_name][outcome_label] = boot_df

            rank_df = bootstrap_rank_stability(
                method, X_arr, Y_arr, factor_names,
                n_bootstrap=args.n_bootstrap, seed=args.seed,
                method_kwargs=kw,
            )
            rank_stab_results[method_name][outcome_label] = rank_df

            conv_df = convergence_curve(
                method, X_arr, Y_arr, factor_names,
                sizes=args.convergence_sizes,
                n_replicates=args.n_convergence_replicates,
                seed=args.seed, method_kwargs=kw,
            )
            convergence_results[method_name][outcome_label] = conv_df

    # --- Cross-method rank-correlation (per outcome) ---
    cm_corr_records = []
    for outcome_label, *_ in full_outcome_list:
        per_method = {
            m: method_results[m][outcome_label] for m in args.methods
        }
        rho = cross_method_rank_corr(per_method)
        for mi in rho.index:
            for mj in rho.columns:
                cm_corr_records.append({
                    "outcome": outcome_label,
                    "method_i": mi,
                    "method_j": mj,
                    "rho": float(rho.loc[mi, mj]),
                })
    cm_corr_df = pd.DataFrame(cm_corr_records)

    # --- Cross-outcome rank-correlation (per method) ---
    co_corr_records = []
    for method_name in args.methods:
        per_outcome = method_results[method_name]
        rho = cross_outcome_rank_corr(per_outcome, method=method_name)
        for oi in rho.index:
            for oj in rho.columns:
                co_corr_records.append({
                    "method": method_name,
                    "outcome_i": oi,
                    "outcome_j": oj,
                    "rho": float(rho.loc[oi, oj]),
                })
    co_corr_df = pd.DataFrame(co_corr_records)

    # --- Selection criterion (per outcome × method) ---
    selection_log: Dict[str, Dict[str, dict]] = {}
    for outcome_label, *_ in full_outcome_list:
        outcome_decisions: Dict[str, dict] = {}
        for method_name in args.methods:
            comparators = [m for m in args.methods if m != method_name]
            cross_rho_for_method = pd.Series({
                comp: float(
                    cm_corr_df.query(
                        "outcome == @outcome_label and method_i == @method_name "
                        "and method_j == @comp"
                    )["rho"].iloc[0]
                )
                for comp in comparators
                if not cm_corr_df.query(
                    "outcome == @outcome_label and method_i == @method_name "
                    "and method_j == @comp"
                ).empty
            })
            decision = apply_method_selection_criterion(
                bootstrap_df=bootstrap_results[method_name][outcome_label],
                rank_stability_df=rank_stab_results[method_name][outcome_label],
                cross_method_rho=cross_rho_for_method,
                rank_spearman_threshold=args.rank_spearman_threshold,
                cross_method_threshold=args.cross_method_threshold,
            )
            outcome_decisions[method_name] = decision

        # Anchor selection: prefer Delta if it passes, else PAWN, else RBD-FAST.
        anchor: Optional[str] = None
        for preferred in ("delta", "pawn", "rbd_fast"):
            if preferred in outcome_decisions and outcome_decisions[preferred]["passes"]:
                anchor = preferred
                break
        outcome_decisions["_anchor"] = anchor or "none"
        selection_log[outcome_label] = outcome_decisions

    # --- Persist long-form result tables ---
    def _long_form(
        per_outcome: Dict[str, pd.DataFrame],
        *,
        headline_col: Optional[str] = None,
    ) -> pd.DataFrame:
        """Stack per-outcome DataFrames into long form.

        If ``headline_col`` is supplied, a ``headline_index`` column is
        added carrying the value of that column (used for the per-method
        indices DataFrames so consumers don't have to look up
        :data:`HEADLINE_INDEX`).
        """
        rows = []
        for outcome_label, df in per_outcome.items():
            for factor in df.index:
                row = {"outcome": outcome_label, "factor": str(factor)}
                row.update({c: float(df.loc[factor, c]) for c in df.columns
                            if pd.api.types.is_numeric_dtype(df[c])})
                if headline_col is not None and headline_col in df.columns:
                    row["headline_index"] = float(df.loc[factor, headline_col])
                rows.append(row)
        return pd.DataFrame(rows)

    for method_name in args.methods:
        headline = HEADLINE_INDEX[method_name]
        idx_long = _long_form(method_results[method_name], headline_col=headline)
        boot_long = _long_form(bootstrap_results[method_name])
        rank_long = _long_form(rank_stab_results[method_name])
        # Stash the rank-Spearman summary attrs alongside the rank table
        # via a small companion frame so the diagnostic layer can recover
        # the threshold-checking value from disk.
        rank_summary = pd.DataFrame([
            {
                "outcome": oc,
                "rank_spearman_median": rank_stab_results[method_name][oc].attrs.get(
                    "rank_spearman_median"
                ),
                "rank_spearman_iqr_lo": rank_stab_results[method_name][oc].attrs.get(
                    "rank_spearman_iqr_lo"
                ),
                "rank_spearman_iqr_hi": rank_stab_results[method_name][oc].attrs.get(
                    "rank_spearman_iqr_hi"
                ),
            }
            for oc in rank_stab_results[method_name]
        ])
        conv_long = pd.concat(
            [df.assign(outcome=oc, method=method_name)
             for oc, df in convergence_results[method_name].items()],
            ignore_index=True,
        )
        idx_long.to_parquet(results_dir / f"indices_{method_name}.parquet")
        boot_long.to_parquet(results_dir / f"bootstrap_{method_name}.parquet")
        rank_long.to_parquet(results_dir / f"rank_stability_{method_name}.parquet")
        rank_summary.to_parquet(
            results_dir / f"rank_stability_summary_{method_name}.parquet"
        )
        conv_long.to_parquet(results_dir / f"convergence_{method_name}.parquet")

    cm_corr_df.to_parquet(results_dir / "cross_method_rank_corr.parquet")
    co_corr_df.to_parquet(results_dir / "cross_outcome_rank_corr.parquet")

    selection_log_payload = {
        "_meta": {
            "rank_spearman_threshold": args.rank_spearman_threshold,
            "cross_method_threshold": args.cross_method_threshold,
            "n_realizations_per_outcome": realization_counts,
        },
        **selection_log,
    }
    (results_dir / "selection_log.json").write_text(
        json.dumps(selection_log_payload, indent=2, default=str)
    )

    print(f"[08] computed all indices in {time.time() - t0:.1f}s")
    print(f"[08] selection anchors per outcome: " +
          ", ".join(f"{oc}={selection_log[oc]['_anchor']}"
                    for oc, *_ in full_outcome_list))
    print(f"[08] done. results -> {results_dir}")


if __name__ == "__main__":
    main()
