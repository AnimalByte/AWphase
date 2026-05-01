#!/usr/bin/env python3
import argparse
import csv
import json
import traceback
from pathlib import Path

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
    HAS_SGK = True
except ImportError:
    from sklearn.model_selection import GroupKFold
    HAS_SGK = False

from xgboost import XGBClassifier


EXCLUDE_COLUMNS = {
    "label",
    "source_window",
    "block_a",
    "block_b",
    "true_relation",
    "block_a_orientation_sites",
    "block_b_orientation_sites",
    "block_a_orientation_margin",
    "block_b_orientation_margin",
    "pred_relation",
    "bridge_pred_relation",
    "block_a_start",
    "block_a_end",
    "block_b_start",
    "block_b_end",
}


def json_safe(x):
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [json_safe(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.ndarray):
        return json_safe(x.tolist())
    if isinstance(x, pd.Series):
        return json_safe(x.tolist())
    if isinstance(x, pd.Index):
        return json_safe(list(x))
    try:
        if pd.isna(x) and not isinstance(x, (str, bytes)):
            return None
    except Exception:
        pass
    return x


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(json_safe(obj), indent=2) + "\n")


def pick_features(df, label_col, group_col):
    excluded = set(EXCLUDE_COLUMNS) | {label_col, group_col}
    features = []
    constant = []
    skipped = []

    for c in df.columns:
        if c in excluded:
            continue

        s = pd.to_numeric(df[c], errors="coerce")
        if s.notna().sum() == 0:
            skipped.append(c)
            continue

        s = s.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        if s.nunique(dropna=True) <= 1:
            constant.append(c)
        else:
            features.append(c)

    if not features:
        raise SystemExit("No usable numeric non-leakage features found.")

    return features, constant, skipped


def build_x(df, features):
    return pd.DataFrame(
        {
            f: pd.to_numeric(df[f], errors="coerce")
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            for f in features
        },
        index=df.index,
    )


def make_model(seed, threads, n_estimators, max_depth, learning_rate, scale_pos_weight):
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
        n_jobs=threads,
        tree_method="hist",
        scale_pos_weight=scale_pos_weight,
    )


def main_inner(args):
    training = Path(args.training_tsv)
    if not training.exists() or training.stat().st_size == 0:
        raise SystemExit(f"Missing or empty training TSV: {training}")

    print(f"[load] {training}")
    df = pd.read_csv(training, sep="\t", low_memory=False)
    print(f"[load] rows={len(df)} cols={len(df.columns)}")

    if len(df) == 0:
        raise SystemExit("Training TSV has zero rows.")

    if args.label_col not in df.columns:
        raise SystemExit(f"Missing label column: {args.label_col}")
    if args.group_col not in df.columns:
        raise SystemExit(f"Missing group column: {args.group_col}")

    y = pd.to_numeric(df[args.label_col], errors="coerce").fillna(0).astype(int).values
    groups = df[args.group_col].astype(str).values

    pos = int((y == 1).sum())
    neg = int((y == 0).sum())

    print(f"[target] positives={pos} negatives={neg} positive_rate={pos / len(df) if len(df) else None}")

    if pos == 0 or neg == 0:
        raise SystemExit(f"Need both classes. positives={pos}, negatives={neg}")

    features, constant, skipped = pick_features(df, args.label_col, args.group_col)
    print(f"[features] usable={len(features)} constant={len(constant)} skipped_non_numeric={len(skipped)}")
    print("[features] selected:")
    for f in features:
        print(f"  {f}")

    feature_report = {
        "training_tsv": str(training),
        "rows": int(len(df)),
        "columns": list(df.columns),
        "features": features,
        "constant_features": constant,
        "skipped_non_numeric_columns": skipped,
        "excluded_columns": sorted(EXCLUDE_COLUMNS),
        "positives": pos,
        "negatives": neg,
    }
    write_json(args.out_feature_report_json, feature_report)

    X = build_x(df, features)

    unique_groups = sorted(set(groups))
    n_splits = min(args.n_splits, len(unique_groups))

    print(f"[cv] unique_groups={len(unique_groups)} n_splits={n_splits}")

    if n_splits < 2:
        raise SystemExit(f"Need at least 2 groups. Found {len(unique_groups)}")

    scale_pos_weight = neg / pos if pos else 1.0
    probs = np.zeros(len(df), dtype=float)
    fold_rows = []

    if HAS_SGK:
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=args.seed)
        split_iter = splitter.split(X, y, groups)
        splitter_name = "StratifiedGroupKFold"
    else:
        splitter = GroupKFold(n_splits=n_splits)
        split_iter = splitter.split(X, y, groups)
        splitter_name = "GroupKFold"

    for fold, (tr, te) in enumerate(split_iter, start=1):
        print(f"[fold {fold}] train={len(tr)} test={len(te)} test_pos={int(y[te].sum())} test_neg={int((y[te] == 0).sum())}")

        model = make_model(
            seed=args.seed + fold,
            threads=args.threads,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            learning_rate=args.learning_rate,
            scale_pos_weight=scale_pos_weight,
        )

        model.fit(X.iloc[tr], y[tr])
        p = model.predict_proba(X.iloc[te])[:, 1]
        probs[te] = p

        yhat = (p >= 0.5).astype(int)
        pr, rc, f1, _ = precision_recall_fscore_support(
            y[te], yhat, average="binary", zero_division=0
        )

        row = {
            "fold": fold,
            "test_n": int(len(te)),
            "test_pos": int(y[te].sum()),
            "test_neg": int((y[te] == 0).sum()),
            "accuracy_0.5": float(accuracy_score(y[te], yhat)),
            "precision_0.5": float(pr),
            "recall_0.5": float(rc),
            "f1_0.5": float(f1),
            "test_groups": sorted(set(groups[te])),
        }

        if len(set(y[te])) == 2:
            row["roc_auc"] = float(roc_auc_score(y[te], p))
            row["average_precision"] = float(average_precision_score(y[te], p))
        else:
            row["roc_auc"] = None
            row["average_precision"] = None

        fold_rows.append(row)

    yhat_all = (probs >= 0.5).astype(int)
    pr, rc, f1, _ = precision_recall_fscore_support(
        y, yhat_all, average="binary", zero_division=0
    )

    metrics = {
        "training_tsv": str(training),
        "rows": int(len(df)),
        "positives": pos,
        "negatives": neg,
        "positive_rate": float(pos / len(df)),
        "label_col": args.label_col,
        "group_col": args.group_col,
        "splitter": splitter_name,
        "n_splits": n_splits,
        "selected_feature_count": len(features),
        "selected_features": features,
        "constant_features": constant,
        "skipped_non_numeric_columns": skipped,
        "excluded_columns": sorted(EXCLUDE_COLUMNS),
        "cv_roc_auc": float(roc_auc_score(y, probs)),
        "cv_average_precision": float(average_precision_score(y, probs)),
        "cv_accuracy_0.5": float(accuracy_score(y, yhat_all)),
        "cv_precision_0.5": float(pr),
        "cv_recall_0.5": float(rc),
        "cv_f1_0.5": float(f1),
        "folds": fold_rows,
    }

    print("[final] fit full model")
    final_model = make_model(
        seed=args.seed,
        threads=args.threads,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        scale_pos_weight=scale_pos_weight,
    )
    final_model.fit(X, y)

    Path(args.out_model_json).parent.mkdir(parents=True, exist_ok=True)
    final_model.save_model(args.out_model_json)
    print(f"[write] model={args.out_model_json}")

    out_pred = df.copy()
    out_pred["xgb_oof_probability"] = probs
    out_pred["xgb_oof_pred_0.5"] = yhat_all
    Path(args.out_predictions_tsv).parent.mkdir(parents=True, exist_ok=True)
    out_pred.to_csv(args.out_predictions_tsv, sep="\t", index=False)
    print(f"[write] predictions={args.out_predictions_tsv}")

    write_json(args.out_metrics_json, metrics)
    print(f"[write] metrics={args.out_metrics_json}")

    booster = final_model.get_booster()
    gain = booster.get_score(importance_type="gain")
    weight = booster.get_score(importance_type="weight")

    imp_rows = []
    for f in features:
        imp_rows.append({
            "feature": f,
            "gain": float(gain.get(f, 0.0)),
            "weight": float(weight.get(f, 0.0)),
        })

    imp_rows.sort(key=lambda r: r["gain"], reverse=True)

    Path(args.out_feature_importance_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_feature_importance_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["feature", "gain", "weight"], delimiter="\t")
        w.writeheader()
        w.writerows(imp_rows)

    print(f"[write] feature_importance={args.out_feature_importance_tsv}")

    print(json.dumps(json_safe({
        "rows": len(df),
        "positives": pos,
        "negatives": neg,
        "features": len(features),
        "cv_roc_auc": metrics["cv_roc_auc"],
        "cv_average_precision": metrics["cv_average_precision"],
        "predictions": args.out_predictions_tsv,
    }), indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--out-model-json", required=True)
    ap.add_argument("--out-metrics-json", required=True)
    ap.add_argument("--out-feature-importance-tsv", required=True)
    ap.add_argument("--out-predictions-tsv", required=True)
    ap.add_argument("--out-feature-report-json", required=True)
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--group-col", default="source_window")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-estimators", type=int, default=250)
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--learning-rate", type=float, default=0.04)
    args = ap.parse_args()

    try:
        main_inner(args)
    except Exception:
        print()
        print("===== PYTHON TRACEBACK =====")
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
