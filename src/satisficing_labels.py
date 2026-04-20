"""Label + classifier layer for the satisficing framework.

This module is the cheap-to-rerun half of the scenario-discovery
pipeline: given a pre-computed metric bank (see
:mod:`src.satisficing_metrics`) and a manifest of candidate satisficing
definitions, produce binary labels, train a gradient-boosted-tree
classifier per definition on drought characteristics, and emit a
summary table plus per-definition artifacts.

The manifest format is YAML; each entry names a metric column in the
bank, a threshold, a comparison direction, and a human-readable ID.
Changing the manifest and re-running this layer does not require the
Pywr-DRB simulation to re-run — that is the whole point of the
separation.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Manifest parsing
# ---------------------------------------------------------------------------


VALID_DIRECTIONS = {"le", "lt", "ge", "gt"}


@dataclass
class SatisficingDefinition:
    """One candidate satisficing rule parsed from the manifest."""

    definition_id: str
    metric: str
    threshold: float
    direction: str
    description: str = ""

    def apply(self, series: pd.Series) -> pd.Series:
        """Return a binary label (1 = satisficing) for each row in *series*."""
        if self.direction == "le":
            y = series <= self.threshold
        elif self.direction == "lt":
            y = series < self.threshold
        elif self.direction == "ge":
            y = series >= self.threshold
        elif self.direction == "gt":
            y = series > self.threshold
        else:
            raise ValueError(
                f"definition {self.definition_id!r}: invalid direction "
                f"{self.direction!r} (must be one of {sorted(VALID_DIRECTIONS)})"
            )
        return y.astype("Int64")


def load_manifest(path: Path) -> List[SatisficingDefinition]:
    """Load a YAML manifest and return a list of definitions."""
    import yaml

    payload = yaml.safe_load(Path(path).read_text())
    if not isinstance(payload, list):
        raise ValueError(
            f"Manifest {path} must be a YAML list of definition mappings; "
            f"got {type(payload).__name__}"
        )
    defs: List[SatisficingDefinition] = []
    seen = set()
    for i, entry in enumerate(payload):
        if not isinstance(entry, dict):
            raise ValueError(f"Manifest entry {i} is not a mapping: {entry!r}")
        missing = {"definition_id", "metric", "threshold", "direction"} - entry.keys()
        if missing:
            raise ValueError(
                f"Manifest entry {i} missing required keys: {sorted(missing)}"
            )
        if entry["definition_id"] in seen:
            raise ValueError(f"Duplicate definition_id: {entry['definition_id']!r}")
        seen.add(entry["definition_id"])
        direction = str(entry["direction"]).lower()
        if direction not in VALID_DIRECTIONS:
            raise ValueError(
                f"Manifest entry {entry['definition_id']!r}: direction "
                f"{direction!r} not in {sorted(VALID_DIRECTIONS)}"
            )
        defs.append(SatisficingDefinition(
            definition_id=str(entry["definition_id"]),
            metric=str(entry["metric"]),
            threshold=float(entry["threshold"]),
            direction=direction,
            description=str(entry.get("description", "")),
        ))
    return defs


# ---------------------------------------------------------------------------
# Label application
# ---------------------------------------------------------------------------


def apply_labels(
    bank_df: pd.DataFrame,
    manifest: Sequence[SatisficingDefinition],
) -> pd.DataFrame:
    """Apply every manifest definition to the bank.

    Returns a long-form DataFrame with columns
    ``[definition_id, realization_id, metric_value, threshold, direction, y]``.
    Rows where the metric is NaN are dropped (the classifier cannot
    consume missing labels).
    """
    bank = bank_df.reset_index()
    rows: List[pd.DataFrame] = []
    for d in manifest:
        if d.metric not in bank.columns:
            warnings.warn(
                f"[satisficing_labels] definition {d.definition_id!r} "
                f"references missing metric {d.metric!r}; skipping",
                stacklevel=2,
            )
            continue
        metric_series = pd.to_numeric(bank[d.metric], errors="coerce")
        labels = d.apply(metric_series)
        chunk = pd.DataFrame({
            "definition_id": d.definition_id,
            "realization_id": bank["realization_id"].astype(str).values,
            "metric": d.metric,
            "metric_value": metric_series.values,
            "threshold": d.threshold,
            "direction": d.direction,
            "y": labels.values,
        })
        rows.append(chunk.dropna(subset=["metric_value", "y"]))
    if not rows:
        return pd.DataFrame(columns=[
            "definition_id", "realization_id", "metric",
            "metric_value", "threshold", "direction", "y",
        ])
    out = pd.concat(rows, ignore_index=True)
    out["y"] = out["y"].astype(int)
    return out


# ---------------------------------------------------------------------------
# Classifier training
# ---------------------------------------------------------------------------


def _fit_classifier_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    make_model,
    importance_fn,
    seed: int = 42,
    n_splits_primary: int = 5,
    min_class_for_primary: int = 15,
    n_splits_fallback: int = 3,
) -> Dict:
    """Shared stratified-K-fold CV + final-fit loop.

    Args:
        make_model: ``make_model(seed) -> fresh sklearn estimator``. Called
            once per fold and once for the final model.
        importance_fn: ``importance_fn(fitted_model, feature_names) -> dict``.
            Returns a per-feature importance score (e.g., ``.feature_importances_``
            for trees, absolute scaled coefficients for logistic regression).

    Falls back from 5-fold to 3-fold when the minority class has fewer than
    ``min_class_for_primary`` members. Returns ``status="degenerate_labels"``
    with NaN scores when either class is empty.
    """
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score, brier_score_loss

    y = np.asarray(y, dtype=int)
    classes, counts = np.unique(y, return_counts=True)
    if len(classes) < 2:
        return {
            "status": "degenerate_labels",
            "class_counts": {int(c): int(n) for c, n in zip(classes, counts)},
            "n_splits": 0,
            "auc_mean": float("nan"),
            "auc_std": float("nan"),
            "brier_mean": float("nan"),
            "brier_std": float("nan"),
            "feature_importance": {},
            "cv_predictions": None,
            "model": None,
        }

    min_class = int(counts.min())
    X_arr = X.values
    feature_names = list(X.columns)

    # With a single-member minority class, stratified CV can't put a
    # non-empty minority in both train and test folds. Train a final
    # model on all data and report no CV metrics.
    if min_class < 2:
        final_model = make_model(seed)
        final_model.fit(X_arr, y)
        return {
            "status": "cv_infeasible_single_member",
            "class_counts": {int(c): int(n) for c, n in zip(classes, counts)},
            "n_splits": 0,
            "auc_mean": float("nan"),
            "auc_std": float("nan"),
            "brier_mean": float("nan"),
            "brier_std": float("nan"),
            "feature_importance": importance_fn(final_model, feature_names),
            "cv_predictions": None,
            "model": final_model,
        }

    n_splits = n_splits_primary if min_class >= min_class_for_primary else n_splits_fallback
    if min_class < n_splits:
        n_splits = max(2, min_class)

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

    # Per-fold predictions so we can score AUC + Brier consistently across
    # the same fold split on pooled predictions and per-fold AUC std.
    cv_preds = np.full(len(y), np.nan, dtype=float)
    cv_labels = np.full(len(y), -1, dtype=int)
    fold_aucs: List[float] = []
    for train_idx, test_idx in skf.split(X_arr, y):
        model_fold = make_model(seed)
        model_fold.fit(X_arr[train_idx], y[train_idx])
        probs = model_fold.predict_proba(X_arr[test_idx])[:, 1]
        cv_preds[test_idx] = probs
        cv_labels[test_idx] = y[test_idx]
        if len(np.unique(y[test_idx])) > 1:
            fold_aucs.append(float(roc_auc_score(y[test_idx], probs)))

    auc = roc_auc_score(cv_labels, cv_preds) if len(np.unique(cv_labels)) > 1 else float("nan")
    brier = brier_score_loss(cv_labels, cv_preds)

    # Final model on all data for feature importance and downstream plotting
    final_model = make_model(seed)
    final_model.fit(X_arr, y)

    status = "ok" if n_splits == n_splits_primary else "cv_reduced_folds"
    return {
        "status": status,
        "class_counts": {int(c): int(n) for c, n in zip(classes, counts)},
        "n_splits": int(n_splits),
        "auc_mean": float(np.nanmean(fold_aucs)) if fold_aucs else float(auc),
        "auc_std": float(np.nanstd(fold_aucs, ddof=1)) if len(fold_aucs) > 1 else float("nan"),
        "brier_mean": float(brier),
        "brier_std": float("nan"),
        "feature_importance": importance_fn(final_model, feature_names),
        "cv_predictions": pd.DataFrame({
            "y_true": cv_labels,
            "y_prob": cv_preds,
        }),
        "model": final_model,
    }


def fit_gbt_classifier(
    X: pd.DataFrame,
    y: np.ndarray,
    seed: int = 42,
    n_splits_primary: int = 5,
    min_class_for_primary: int = 15,
    n_splits_fallback: int = 3,
) -> Dict:
    """Fit a gradient-boosted-tree classifier with stratified K-fold CV.

    Falls back from 5-fold to 3-fold when the minority class has fewer than
    *min_class_for_primary* members. Returns ``status="degenerate_labels"``
    with NaN scores when either class is empty.
    """
    from sklearn.ensemble import GradientBoostingClassifier

    def _make(s):
        return GradientBoostingClassifier(random_state=s)

    def _importance(model, names):
        return dict(zip(
            names,
            [float(v) for v in model.feature_importances_],
        ))

    return _fit_classifier_cv(
        X, y, _make, _importance,
        seed=seed,
        n_splits_primary=n_splits_primary,
        min_class_for_primary=min_class_for_primary,
        n_splits_fallback=n_splits_fallback,
    )


def fit_logreg_classifier(
    X: pd.DataFrame,
    y: np.ndarray,
    seed: int = 42,
    poly_degree: int = 2,
    n_splits_primary: int = 5,
    min_class_for_primary: int = 15,
    n_splits_fallback: int = 3,
) -> Dict:
    """Fit a logistic-regression classifier with polynomial features.

    The pipeline is ``PolynomialFeatures(degree=poly_degree) → StandardScaler
    → LogisticRegression`` (the scaler/bias structure matches the standard
    scenario-discovery convention; quadratic features reproduce the curved
    boundaries common in the SD literature while preserving the smooth,
    calibrated probabilities of a GLM). Feature importance is reported as
    the absolute value of the standardized coefficients on the original
    (degree-1) features; higher-order terms contribute to the fitted
    probabilities but are not surfaced individually in the summary table.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler, PolynomialFeatures
    from sklearn.pipeline import Pipeline

    def _make(s):
        steps = []
        if poly_degree > 1:
            steps.append((
                "poly",
                PolynomialFeatures(degree=poly_degree, include_bias=False),
            ))
        steps.append(("scaler", StandardScaler()))
        steps.append((
            "logreg",
            LogisticRegression(random_state=s, max_iter=1000),
        ))
        return Pipeline(steps)

    def _importance(model, names):
        coef = np.asarray(model.named_steps["logreg"].coef_[0], dtype=float)
        n = len(names)
        return dict(zip(names, [float(abs(v)) for v in coef[:n]]))

    return _fit_classifier_cv(
        X, y, _make, _importance,
        seed=seed,
        n_splits_primary=n_splits_primary,
        min_class_for_primary=min_class_for_primary,
        n_splits_fallback=n_splits_fallback,
    )


_MODEL_FITTERS = {
    "gbt": fit_gbt_classifier,
    "logreg": fit_logreg_classifier,
}


def fit_classifier(
    X: pd.DataFrame,
    y: np.ndarray,
    model_name: str = "gbt",
    **kwargs,
) -> Dict:
    """Dispatch to a named classifier fitter (``gbt`` or ``logreg``)."""
    if model_name not in _MODEL_FITTERS:
        raise ValueError(
            f"Unknown classifier model_name={model_name!r}; "
            f"expected one of {list(_MODEL_FITTERS)}"
        )
    return _MODEL_FITTERS[model_name](X, y, **kwargs)


# ---------------------------------------------------------------------------
# Manifest sweep
# ---------------------------------------------------------------------------


def sweep_manifest(
    bank_df: pd.DataFrame,
    chars_df: pd.DataFrame,
    manifest: Sequence[SatisficingDefinition],
    feature_cols: Sequence[str] = ("mean_duration", "mean_avg_severity", "peak_severity_month"),
    output_dir: Optional[Path] = None,
    seed: int = 42,
    model_name: str = "gbt",
) -> pd.DataFrame:
    """Run every manifest definition through label → classifier.

    Args:
        bank_df: Metric bank from :func:`src.satisficing_metrics.compute_metric_bank`.
        chars_df: Drought characteristics indexed by realization_id (the
            *features* passed to the classifier).
        manifest: List of :class:`SatisficingDefinition`.
        feature_cols: Drought-characteristic columns in *chars_df* used as
            classifier features.
        output_dir: If given, per-definition artifacts (``model.joblib``,
            ``cv_predictions.csv``) are written under
            ``output_dir/classifiers/{definition_id}/``.
        seed: Shared random seed for classifier and CV folds.
        model_name: Which classifier to dispatch — ``"gbt"`` (default,
            gradient-boosted trees) or ``"logreg"`` (logistic regression
            with degree-2 polynomial features). See :func:`fit_classifier`.

    Returns:
        Summary DataFrame with one row per definition.
    """
    import joblib

    labels_long = apply_labels(bank_df, manifest)
    summary_rows: List[Dict] = []

    for d in manifest:
        sub = labels_long[labels_long["definition_id"] == d.definition_id]
        if sub.empty:
            summary_rows.append({
                "definition_id": d.definition_id,
                "metric": d.metric,
                "threshold": d.threshold,
                "direction": d.direction,
                "description": d.description,
                "status": "no_labels",
                "n_pos": 0, "n_neg": 0,
                "auc_mean": float("nan"), "auc_std": float("nan"),
                "brier_mean": float("nan"), "n_splits": 0,
                **{f"feat_imp_{c}": float("nan") for c in feature_cols},
            })
            continue

        real_ids = sub["realization_id"].astype(str).values
        missing = [r for r in real_ids if r not in chars_df.index]
        if missing:
            warnings.warn(
                f"[satisficing_labels] definition {d.definition_id!r}: "
                f"{len(missing)} realization_ids not in drought characteristics; "
                f"dropping",
                stacklevel=2,
            )
        keep_mask = np.array([r in chars_df.index for r in real_ids])
        real_ids = real_ids[keep_mask]
        sub = sub.iloc[keep_mask]

        missing_feats = [c for c in feature_cols if c not in chars_df.columns]
        if missing_feats:
            raise KeyError(
                f"drought-characteristics DataFrame missing expected "
                f"feature columns {missing_feats}"
            )
        X = chars_df.loc[real_ids, list(feature_cols)].astype(float)
        y = sub["y"].values.astype(int)

        result = fit_classifier(X, y, model_name=model_name, seed=seed)
        n_pos = int(result["class_counts"].get(1, 0))
        n_neg = int(result["class_counts"].get(0, 0))

        row = {
            "definition_id": d.definition_id,
            "metric": d.metric,
            "threshold": d.threshold,
            "direction": d.direction,
            "description": d.description,
            "status": result["status"],
            "n_pos": n_pos,
            "n_neg": n_neg,
            "auc_mean": result["auc_mean"],
            "auc_std": result["auc_std"],
            "brier_mean": result["brier_mean"],
            "n_splits": result["n_splits"],
        }
        for c in feature_cols:
            row[f"feat_imp_{c}"] = result["feature_importance"].get(c, float("nan"))
        summary_rows.append(row)

        if output_dir is not None and result["model"] is not None:
            out_def = Path(output_dir) / "classifiers" / d.definition_id
            out_def.mkdir(parents=True, exist_ok=True)
            joblib.dump(result["model"], out_def / "model.joblib")
            if result["cv_predictions"] is not None:
                result["cv_predictions"].to_csv(out_def / "cv_predictions.csv",
                                                index=False)
            (out_def / "metadata.json").write_text(json.dumps({
                "definition_id": d.definition_id,
                "metric": d.metric,
                "threshold": d.threshold,
                "direction": d.direction,
                "status": result["status"],
                "class_counts": result["class_counts"],
                "feature_importance": result["feature_importance"],
            }, indent=2))

    return pd.DataFrame(summary_rows)
