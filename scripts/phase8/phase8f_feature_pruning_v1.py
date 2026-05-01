#!/usr/bin/env python3
"""
Phase8F feature pruning for AWPhase candidate gating.

Does two passes:
  1. Greedy family/group pruning
  2. Greedy individual feature pruning after family pruning

Uses grouped cross-validation by physical window so rows from the same physical
window do not leak between train/test folds.
"""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    precision_recall_fscore_support,
    roc_auc_score,
)

try:
    from sklearn.model_selection import StratifiedGroupKFold
    HAS_STRATIFIED_GROUP_KFOLD = True
except ImportError:
    from sklearn.model_selection import GroupKFold
    HAS_STRATIFIED_GROUP_KFOLD = False

from xgboost import XGBClassifier


ID_LIKE_COLUMNS = {
    "source_window",
    "physical_group",
    "chrom",
    "pos",
    "position",
    "variant_id",
    "block_id",
    "selected_block_id",
    "pred_state",
    "truth_state",
    "block_orientation",
    "oriented_pred_state",
    "orientation_same_count",
    "orientation_opposite_count",
    "orientation_margin",
}

TARGET_CANDIDATES = [
    "target",
    "y",
    "is_correct",
    "correct",
    "accepted",
    "label",
    "truth_correct",
    "correct_after_block_orientation",
]


def physical_group_from_source_window(x: str) -> str:
    x = str(x)
    for suffix in [
        "_illumina30x_phase7a",
        "_illumina35x_phase7a",
        "_phase7a_general_runner",
        "_phase7a",
        "_phase8f",
        "_phase8e",
    ]:
        if x.endswith(suffix):
            return x[: -len(suffix)]
    return x


def to_binary_series(s: pd.Series) -> Optional[pd.Series]:
    """Return 0/1 integer series if the column is binary-like, else None."""
    if s.dtype == bool:
        return s.astype(int)

    mapped = s.map(
        lambda x: {
            True: 1,
            False: 0,
            "true": 1,
            "false": 0,
            "TRUE": 1,
            "FALSE": 0,
            "yes": 1,
            "no": 0,
            "YES": 1,
            "NO": 0,
        }.get(x, x)
    )
    numeric = pd.to_numeric(mapped, errors="coerce")
    non_null = numeric.dropna()
    if non_null.empty:
        return None

    values = (
        set(non_null.astype(int).unique())
        if np.all(np.isclose(non_null, non_null.astype(int)))
        else set(non_null.unique())
    )
    if values.issubset({0, 1}) and len(values) >= 1:
        return numeric.fillna(0).astype(int)
    return None


def infer_label_col(df: pd.DataFrame, explicit_label_col: Optional[str] = None) -> str:
    if explicit_label_col:
        if explicit_label_col not in df.columns:
            raise SystemExit(f"Requested --label-col '{explicit_label_col}' is not in the TSV.")
        binary = to_binary_series(df[explicit_label_col])
        if binary is None:
            raise SystemExit(f"Requested --label-col '{explicit_label_col}' is not binary-like.")
        df[explicit_label_col] = binary
        return explicit_label_col

    for c in TARGET_CANDIDATES:
        if c in df.columns:
            binary = to_binary_series(df[c])
            if binary is not None:
                df[c] = binary
                return c

    raise SystemExit(
        "Could not infer binary target column. Use --label-col explicitly. "
        f"Checked: {', '.join(TARGET_CANDIDATES)}"
    )


def infer_group_col(df: pd.DataFrame, explicit_group_col: Optional[str] = None) -> str:
    if explicit_group_col:
        if explicit_group_col not in df.columns:
            raise SystemExit(f"Requested --group-col '{explicit_group_col}' is not in the TSV.")
        return explicit_group_col

    if "physical_group" in df.columns:
        return "physical_group"

    if "source_window" in df.columns:
        df["physical_group"] = df["source_window"].map(physical_group_from_source_window)
        return "physical_group"

    raise SystemExit("Could not infer group column. Expected source_window or physical_group, or pass --group-col.")


def numeric_feature_candidates(df: pd.DataFrame, label_col: str, group_col: str) -> Tuple[List[str], List[str], List[str]]:
    usable: List[str] = []
    constant: List[str] = []
    skipped: List[str] = []

    excluded = set(ID_LIKE_COLUMNS) | {label_col, group_col}

    for c in df.columns:
        if c in excluded:
            continue

        s = pd.to_numeric(df[c], errors="coerce")
        non_null = int(s.notna().sum())
        if non_null == 0:
            skipped.append(c)
            continue

        s = s.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        if s.nunique(dropna=True) <= 1:
            constant.append(c)
        else:
            usable.append(c)

    if not usable:
        raise SystemExit("No usable numeric features found after excluding IDs, target, group, and constants.")

    return usable, constant, skipped


def build_X(df: pd.DataFrame, features: List[str]) -> pd.DataFrame:
    if not features:
        raise ValueError("Cannot build feature matrix with zero features.")

    return pd.DataFrame(
        {
            f: pd.to_numeric(df[f], errors="coerce")
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            for f in features
        },
        index=df.index,
    )


def metric_key(features: List[str]) -> str:
    joined = "\n".join(sorted(features))
    return hashlib.md5(joined.encode()).hexdigest()


def safe_delta(a, b):
    if a is None or b is None:
        return None
    return float(a - b)


def jsonable_metric_dict(d: dict) -> dict:
    out = {}
    for k, v in d.items():
        if k == "_oof_probability":
            continue
        if isinstance(v, np.ndarray):
            continue
        if isinstance(v, (np.integer,)):
            out[k] = int(v)
        elif isinstance(v, (np.floating,)):
            out[k] = float(v)
        else:
            out[k] = v
    return out


def make_classifier(
    *,
    n_estimators: int,
    max_depth: int,
    learning_rate: float,
    seed: int,
    n_threads: int,
    scale_pos_weight: float,
) -> XGBClassifier:
    return XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=2,
        reg_lambda=2.0,
        reg_alpha=0.25,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=seed,
        n_jobs=n_threads,
        tree_method="hist",
        scale_pos_weight=scale_pos_weight,
    )


def evaluate_feature_set(
    *,
    df: pd.DataFrame,
    features: List[str],
    label_col: str,
    group_col: str,
    n_splits: int,
    n_threads: int,
    seed: int,
    n_estimators: int,
    max_depth: int,
    learning_rate: float,
    cache: Dict[str, dict],
) -> dict:
    features = list(dict.fromkeys(features))
    if not features:
        raise ValueError("evaluate_feature_set received zero features.")

    key = metric_key(features)
    if key in cache:
        return cache[key]

    X = build_X(df, features)
    y = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int).values
    groups = df[group_col].astype(str).values

    unique_groups = sorted(set(groups))
    splits = min(n_splits, len(unique_groups))
    if splits < 2:
        raise SystemExit(f"Need at least 2 groups for grouped CV; found {len(unique_groups)}.")

    pos = int((y == 1).sum())
    neg = int((y == 0).sum())
    if pos == 0 or neg == 0:
        raise SystemExit(f"Target column '{label_col}' must contain both classes. Found positives={pos}, negatives={neg}.")

    scale_pos_weight = neg / pos if pos else 1.0
    probs = np.zeros(len(df), dtype=float)
    fold_rows = []

    if HAS_STRATIFIED_GROUP_KFOLD:
        splitter = StratifiedGroupKFold(n_splits=splits, shuffle=True, random_state=seed)
        split_iter = splitter.split(X, y, groups)
    else:
        splitter = GroupKFold(n_splits=splits)
        split_iter = splitter.split(X, y, groups)

    for fold, (tr, te) in enumerate(split_iter, start=1):
        model = make_classifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            seed=seed + fold,
            n_threads=n_threads,
            scale_pos_weight=scale_pos_weight,
        )
        model.fit(X.iloc[tr], y[tr])
        p = model.predict_proba(X.iloc[te])[:, 1]
        probs[te] = p

        yhat = (p >= 0.5).astype(int)
        pr, rc, f1, _ = precision_recall_fscore_support(y[te], yhat, average="binary", zero_division=0)

        fold_row = {
            "fold": int(fold),
            "test_groups": [str(g) for g in sorted(set(groups[te]))],
            "test_n": int(len(te)),
            "test_pos": int(y[te].sum()),
            "test_neg": int((y[te] == 0).sum()),
            "accuracy_0.5": float(accuracy_score(y[te], yhat)),
            "precision_0.5": float(pr),
            "recall_0.5": float(rc),
            "f1_0.5": float(f1),
        }

        if len(set(y[te])) == 2:
            fold_row["roc_auc"] = float(roc_auc_score(y[te], p))
            fold_row["average_precision"] = float(average_precision_score(y[te], p))
        else:
            fold_row["roc_auc"] = None
            fold_row["average_precision"] = None

        fold_rows.append(fold_row)

    yhat_all = (probs >= 0.5).astype(int)
    pr, rc, f1, _ = precision_recall_fscore_support(y, yhat_all, average="binary", zero_division=0)

    out = {
        "n_features": int(len(features)),
        "features": features,
        "rows": int(len(df)),
        "positives": pos,
        "negatives": neg,
        "positive_rate": float(pos / len(df)) if len(df) else 0.0,
        "cv_roc_auc": float(roc_auc_score(y, probs)) if len(set(y)) == 2 else None,
        "cv_average_precision": float(average_precision_score(y, probs)),
        "cv_accuracy_0.5": float(accuracy_score(y, yhat_all)),
        "cv_precision_0.5": float(pr),
        "cv_recall_0.5": float(rc),
        "cv_f1_0.5": float(f1),
        "folds": fold_rows,
        "_oof_probability": probs,
    }
    cache[key] = out
    return out


def default_feature_families(features: List[str]) -> Dict[str, List[str]]:
    fam: Dict[str, List[str]] = {
        "rule": [],
        "panel_raw": [],
        "panel_ratio_margin": [],
        "panel_hap_anchor": [],
        "phase8e_population_prior": [],
        "phase8e_weighted_support": [],
        "phase8e_unweighted_support": [],
        "phase8e_logs": [],
        "phase8e_other": [],
        "phase8f": [],
        "other": [],
    }

    for f in features:
        if f == "candidate_accepted_by_rule":
            fam["rule"].append(f)
        elif f in {"panel_confidence", "panel_margin", "panel_support", "panel_conflict"}:
            fam["panel_raw"].append(f)
        elif f in {"support_conflict_ratio", "best_vs_second_margin"}:
            fam["panel_ratio_margin"].append(f)
        elif f in {"panel_haplotypes", "anchors", "n_block_candidates"}:
            fam["panel_hap_anchor"].append(f)
        elif f.startswith("phase8e_pop_") or f in {"phase8e_pop_effective_n", "log_phase8e_pop_effective_n"}:
            fam["phase8e_population_prior"].append(f)
        elif f.startswith("phase8e_weighted_"):
            fam["phase8e_weighted_support"].append(f)
        elif f.startswith("phase8e_unweighted_"):
            fam["phase8e_unweighted_support"].append(f)
        elif f.startswith("log_phase8e_"):
            fam["phase8e_logs"].append(f)
        elif f.startswith("phase8e_"):
            fam["phase8e_other"].append(f)
        elif f.startswith("phase8f_"):
            fam["phase8f"].append(f)
        else:
            fam["other"].append(f)

    return {k: v for k, v in fam.items() if v}


def score_pass(candidate: dict, baseline: dict, metric: str, tolerance: float, ap_tolerance: float) -> bool:
    base_metric = baseline.get(metric)
    cand_metric = candidate.get(metric)
    base_ap = baseline.get("cv_average_precision")
    cand_ap = candidate.get("cv_average_precision")

    if cand_metric is None or base_metric is None or base_ap is None or cand_ap is None:
        return False

    return cand_metric >= base_metric - tolerance and cand_ap >= base_ap - ap_tolerance


def write_rows_tsv(path: Path, rows: List[dict], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def train_final_oof(
    *,
    df: pd.DataFrame,
    features: List[str],
    label_col: str,
    group_col: str,
    args,
    out_predictions_tsv: Path,
    out_model_json: Path,
    out_metrics_json: Path,
    out_feature_importance_tsv: Path,
) -> dict:
    cache: Dict[str, dict] = {}
    m = evaluate_feature_set(
        df=df,
        features=features,
        label_col=label_col,
        group_col=group_col,
        n_splits=args.n_splits,
        n_threads=args.threads,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        cache=cache,
    )

    X = build_X(df, features)
    y = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int).values
    pos = int((y == 1).sum())
    neg = int((y == 0).sum())
    scale_pos_weight = neg / pos if pos else 1.0

    model = make_classifier(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        seed=args.seed,
        n_threads=args.threads,
        scale_pos_weight=scale_pos_weight,
    )
    model.fit(X, y)

    out_model_json.parent.mkdir(parents=True, exist_ok=True)
    model.save_model(str(out_model_json))

    pred_df = df.copy()
    pred_df["xgb_oof_probability"] = m["_oof_probability"]
    pred_df["xgb_oof_pred_0.5"] = (m["_oof_probability"] >= 0.5).astype(int)
    pred_df.to_csv(out_predictions_tsv, sep="\t", index=False)

    metrics = jsonable_metric_dict(m)
    metrics["selected_features"] = features
    metrics["selected_feature_count"] = len(features)
    out_metrics_json.write_text(json.dumps(metrics, indent=2) + "\n")

    booster = model.get_booster()
    gain = booster.get_score(importance_type="gain")
    weight = booster.get_score(importance_type="weight")
    rows = [
        {
            "feature": f,
            "gain": float(gain.get(f, 0.0)),
            "weight": float(weight.get(f, 0.0)),
        }
        for f in features
    ]
    rows.sort(key=lambda r: float(r["gain"]), reverse=True)
    write_rows_tsv(out_feature_importance_tsv, rows, ["feature", "gain", "weight"])

    return metrics


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label-col", default=None, help="Optional explicit binary target column.")
    ap.add_argument("--group-col", default=None, help="Optional explicit grouped-CV column.")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-estimators", type=int, default=250)
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--learning-rate", type=float, default=0.04)
    ap.add_argument("--metric", choices=["cv_roc_auc", "cv_average_precision"], default="cv_roc_auc")
    ap.add_argument("--group-tolerance", type=float, default=0.001)
    ap.add_argument("--individual-tolerance", type=float, default=0.0005)
    ap.add_argument("--ap-tolerance", type=float, default=0.001)
    ap.add_argument("--max-individual-rounds", type=int, default=100)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.training_tsv, sep="\t")
    label_col = infer_label_col(df, args.label_col)
    group_col = infer_group_col(df, args.group_col)
    features, constant_features, skipped_non_numeric = numeric_feature_candidates(df, label_col, group_col)
    families = default_feature_families(features)

    cache: Dict[str, dict] = {}
    baseline = evaluate_feature_set(
        df=df,
        features=features,
        label_col=label_col,
        group_col=group_col,
        n_splits=args.n_splits,
        n_threads=args.threads,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        cache=cache,
    )

    family_rows = []
    for fam_name, fam_features in sorted(families.items()):
        fam_set = set(fam_features)
        cand_features = [f for f in features if f not in fam_set]
        if not cand_features:
            continue
        m = evaluate_feature_set(
            df=df,
            features=cand_features,
            label_col=label_col,
            group_col=group_col,
            n_splits=args.n_splits,
            n_threads=args.threads,
            seed=args.seed,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            learning_rate=args.learning_rate,
            cache=cache,
        )
        family_rows.append(
            {
                "family_dropped": fam_name,
                "features_in_family": ",".join(fam_features),
                "n_features_after_drop": len(cand_features),
                "cv_roc_auc": m["cv_roc_auc"],
                "cv_average_precision": m["cv_average_precision"],
                "delta_roc_auc": safe_delta(m["cv_roc_auc"], baseline["cv_roc_auc"]),
                "delta_average_precision": safe_delta(m["cv_average_precision"], baseline["cv_average_precision"]),
                "drop_passes_tolerance": score_pass(m, baseline, args.metric, args.group_tolerance, args.ap_tolerance),
            }
        )

    sort_field = "delta_roc_auc" if args.metric == "cv_roc_auc" else "delta_average_precision"
    family_rows.sort(key=lambda r: float(r[sort_field]) if r[sort_field] is not None else -999.0, reverse=True)
    write_rows_tsv(
        outdir / "family_ablation.tsv",
        family_rows,
        [
            "family_dropped",
            "features_in_family",
            "n_features_after_drop",
            "cv_roc_auc",
            "cv_average_precision",
            "delta_roc_auc",
            "delta_average_precision",
            "drop_passes_tolerance",
        ],
    )

    current_features = list(features)
    group_pruning_steps = []
    while True:
        current_base = evaluate_feature_set(
            df=df,
            features=current_features,
            label_col=label_col,
            group_col=group_col,
            n_splits=args.n_splits,
            n_threads=args.threads,
            seed=args.seed,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            learning_rate=args.learning_rate,
            cache=cache,
        )

        candidates = []
        for fam_name, fam_features in sorted(families.items()):
            fam_set = set(fam_features)
            if not any(f in current_features for f in fam_set):
                continue
            cand_features = [f for f in current_features if f not in fam_set]
            if not cand_features:
                continue
            m = evaluate_feature_set(
                df=df,
                features=cand_features,
                label_col=label_col,
                group_col=group_col,
                n_splits=args.n_splits,
                n_threads=args.threads,
                seed=args.seed,
                n_estimators=args.n_estimators,
                max_depth=args.max_depth,
                learning_rate=args.learning_rate,
                cache=cache,
            )
            if score_pass(m, current_base, args.metric, args.group_tolerance, args.ap_tolerance):
                candidates.append((fam_name, fam_features, cand_features, m))

        if not candidates:
            break

        candidates.sort(key=lambda x: (x[3][args.metric], x[3]["cv_average_precision"], -len(x[2])), reverse=True)
        fam_name, fam_features, next_features, next_m = candidates[0]
        group_pruning_steps.append(
            {
                "step": len(group_pruning_steps) + 1,
                "dropped_family": fam_name,
                "dropped_features": ",".join([f for f in fam_features if f in current_features]),
                "n_features_before": len(current_features),
                "n_features_after": len(next_features),
                "before_cv_roc_auc": current_base["cv_roc_auc"],
                "after_cv_roc_auc": next_m["cv_roc_auc"],
                "delta_roc_auc": safe_delta(next_m["cv_roc_auc"], current_base["cv_roc_auc"]),
                "before_cv_average_precision": current_base["cv_average_precision"],
                "after_cv_average_precision": next_m["cv_average_precision"],
                "delta_average_precision": safe_delta(next_m["cv_average_precision"], current_base["cv_average_precision"]),
            }
        )
        current_features = next_features

    write_rows_tsv(
        outdir / "group_pruning_steps.tsv",
        group_pruning_steps,
        [
            "step",
            "dropped_family",
            "dropped_features",
            "n_features_before",
            "n_features_after",
            "before_cv_roc_auc",
            "after_cv_roc_auc",
            "delta_roc_auc",
            "before_cv_average_precision",
            "after_cv_average_precision",
            "delta_average_precision",
        ],
    )

    features_after_group = list(current_features)

    individual_ablation_rows = []
    group_base = evaluate_feature_set(
        df=df,
        features=features_after_group,
        label_col=label_col,
        group_col=group_col,
        n_splits=args.n_splits,
        n_threads=args.threads,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        cache=cache,
    )

    for f in features_after_group:
        cand_features = [x for x in features_after_group if x != f]
        if not cand_features:
            continue
        m = evaluate_feature_set(
            df=df,
            features=cand_features,
            label_col=label_col,
            group_col=group_col,
            n_splits=args.n_splits,
            n_threads=args.threads,
            seed=args.seed,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            learning_rate=args.learning_rate,
            cache=cache,
        )
        individual_ablation_rows.append(
            {
                "feature_dropped": f,
                "n_features_after_drop": len(cand_features),
                "cv_roc_auc": m["cv_roc_auc"],
                "cv_average_precision": m["cv_average_precision"],
                "delta_roc_auc_vs_group_base": safe_delta(m["cv_roc_auc"], group_base["cv_roc_auc"]),
                "delta_average_precision_vs_group_base": safe_delta(
                    m["cv_average_precision"], group_base["cv_average_precision"]
                ),
                "drop_passes_tolerance": score_pass(m, group_base, args.metric, args.individual_tolerance, args.ap_tolerance),
            }
        )

    individual_ablation_rows.sort(
        key=lambda r: float(r["delta_roc_auc_vs_group_base"])
        if r["delta_roc_auc_vs_group_base"] is not None
        else -999.0,
        reverse=True,
    )
    write_rows_tsv(
        outdir / "individual_ablation_after_group.tsv",
        individual_ablation_rows,
        [
            "feature_dropped",
            "n_features_after_drop",
            "cv_roc_auc",
            "cv_average_precision",
            "delta_roc_auc_vs_group_base",
            "delta_average_precision_vs_group_base",
            "drop_passes_tolerance",
        ],
    )

    current_features = list(features_after_group)
    individual_steps = []
    for _ in range(args.max_individual_rounds):
        current_base = evaluate_feature_set(
            df=df,
            features=current_features,
            label_col=label_col,
            group_col=group_col,
            n_splits=args.n_splits,
            n_threads=args.threads,
            seed=args.seed,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            learning_rate=args.learning_rate,
            cache=cache,
        )

        candidates = []
        for f in current_features:
            cand_features = [x for x in current_features if x != f]
            if not cand_features:
                continue
            m = evaluate_feature_set(
                df=df,
                features=cand_features,
                label_col=label_col,
                group_col=group_col,
                n_splits=args.n_splits,
                n_threads=args.threads,
                seed=args.seed,
                n_estimators=args.n_estimators,
                max_depth=args.max_depth,
                learning_rate=args.learning_rate,
                cache=cache,
            )
            if score_pass(m, current_base, args.metric, args.individual_tolerance, args.ap_tolerance):
                candidates.append((f, cand_features, m))

        if not candidates:
            break

        candidates.sort(key=lambda x: (x[2][args.metric], x[2]["cv_average_precision"], -len(x[1])), reverse=True)
        dropped, next_features, next_m = candidates[0]
        individual_steps.append(
            {
                "step": len(individual_steps) + 1,
                "dropped_feature": dropped,
                "n_features_before": len(current_features),
                "n_features_after": len(next_features),
                "before_cv_roc_auc": current_base["cv_roc_auc"],
                "after_cv_roc_auc": next_m["cv_roc_auc"],
                "delta_roc_auc": safe_delta(next_m["cv_roc_auc"], current_base["cv_roc_auc"]),
                "before_cv_average_precision": current_base["cv_average_precision"],
                "after_cv_average_precision": next_m["cv_average_precision"],
                "delta_average_precision": safe_delta(next_m["cv_average_precision"], current_base["cv_average_precision"]),
            }
        )
        current_features = next_features

    write_rows_tsv(
        outdir / "individual_pruning_steps.tsv",
        individual_steps,
        [
            "step",
            "dropped_feature",
            "n_features_before",
            "n_features_after",
            "before_cv_roc_auc",
            "after_cv_roc_auc",
            "delta_roc_auc",
            "before_cv_average_precision",
            "after_cv_average_precision",
            "delta_average_precision",
        ],
    )

    selected_features = list(current_features)
    (outdir / "selected_features.txt").write_text("\n".join(selected_features) + "\n")
    (outdir / "features_initial.txt").write_text("\n".join(features) + "\n")
    (outdir / "features_after_group.txt").write_text("\n".join(features_after_group) + "\n")
    (outdir / "constant_features.txt").write_text("\n".join(constant_features) + "\n")
    (outdir / "skipped_non_numeric_columns.txt").write_text("\n".join(skipped_non_numeric) + "\n")

    final_metrics = train_final_oof(
        df=df,
        features=selected_features,
        label_col=label_col,
        group_col=group_col,
        args=args,
        out_predictions_tsv=outdir / "pruned_oof_predictions.tsv",
        out_model_json=outdir / "pruned_model.json",
        out_metrics_json=outdir / "pruned_metrics.json",
        out_feature_importance_tsv=outdir / "pruned_feature_importance.tsv",
    )

    summary = {
        "training_tsv": args.training_tsv,
        "rows": int(len(df)),
        "label_col": label_col,
        "group_col": group_col,
        "threads": args.threads,
        "n_splits": args.n_splits,
        "splitter": "StratifiedGroupKFold" if HAS_STRATIFIED_GROUP_KFOLD else "GroupKFold",
        "initial_feature_count": len(features),
        "constant_feature_count": len(constant_features),
        "skipped_non_numeric_count": len(skipped_non_numeric),
        "constant_features": constant_features,
        "skipped_non_numeric_columns": skipped_non_numeric,
        "families": families,
        "baseline_without_constants": jsonable_metric_dict(baseline),
        "features_after_group_count": len(features_after_group),
        "selected_feature_count": len(selected_features),
        "selected_features": selected_features,
        "group_pruning_steps": group_pruning_steps,
        "individual_pruning_steps": individual_steps,
        "final_pruned_metrics": final_metrics,
    }
    (outdir / "pruning_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(
        json.dumps(
            {
                "rows": int(len(df)),
                "label_col": label_col,
                "group_col": group_col,
                "splitter": summary["splitter"],
                "initial_features": len(features),
                "constant_features": len(constant_features),
                "skipped_non_numeric_columns": len(skipped_non_numeric),
                "features_after_group": len(features_after_group),
                "selected_features": len(selected_features),
                "baseline_auc": baseline["cv_roc_auc"],
                "baseline_ap": baseline["cv_average_precision"],
                "final_auc": final_metrics["cv_roc_auc"],
                "final_ap": final_metrics["cv_average_precision"],
                "outdir": str(outdir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
