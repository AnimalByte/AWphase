#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import Counter

import numpy as np

try:
    import pandas as pd
    from sklearn.metrics import (
        roc_auc_score,
        average_precision_score,
        accuracy_score,
        precision_recall_fscore_support,
    )
    from sklearn.model_selection import GroupKFold
    from xgboost import XGBClassifier
except Exception as e:
    raise SystemExit(f"Missing dependency or import failure: {e}")

FEATURES = [
    "candidate_accepted_by_rule",
    "panel_confidence",
    "panel_margin",
    "panel_support",
    "panel_conflict",
    "support_conflict_ratio",
    "panel_samples",
    "panel_haplotypes",
    "anchors",
    "best_vs_second_margin",
    "n_block_candidates",
]

def physical_group(source_window: str) -> str:
    """
    Collapse 30x/35x versions of the same physical region into one group.
    This avoids pretending chr20 30x and chr20 35x are independent windows.
    """
    s = source_window
    s = s.replace("_illumina30x_phase7a", "")
    s = s.replace("_illumina35x_phase7a", "")
    s = s.replace("_phase7a_general_runner", "")
    return s

def safe_float(x):
    try:
        if x is None or str(x).strip() == "":
            return 0.0
        return float(x)
    except Exception:
        return 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--out-model-json", required=True)
    ap.add_argument("--out-metrics-json", required=True)
    ap.add_argument("--out-feature-importance-tsv", required=True)
    ap.add_argument("--out-predictions-tsv", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.training_tsv), delimiter="\t"))
    if not rows:
        raise SystemExit("No rows found in training table.")

    missing = [f for f in FEATURES + ["label", "source_window"] if f not in rows[0]]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    y = np.array([int(r["label"]) for r in rows], dtype=int)
    groups = np.array([physical_group(r["source_window"]) for r in rows])
    source_windows = [r["source_window"] for r in rows]

    X = pd.DataFrame(
        [{f: safe_float(r.get(f)) for f in FEATURES} for r in rows],
        columns=FEATURES,
    )

    label_counts = Counter(y)
    group_counts = Counter(groups)
    source_counts = Counter(source_windows)

    if len(label_counts) < 2:
        raise SystemExit(f"Need both positive and negative labels. label_counts={dict(label_counts)}")

    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        raise SystemExit(f"Need at least 2 physical groups. groups={unique_groups}")

    n_splits = min(5, len(unique_groups))

    # Positives are majority here, so this is <1. That is okay.
    pos = int((y == 1).sum())
    neg = int((y == 0).sum())
    scale_pos_weight = neg / pos if pos else 1.0

    params = dict(
        n_estimators=150,
        max_depth=3,
        learning_rate=0.04,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=2,
        reg_lambda=2.0,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=42,
        n_jobs=4,
        tree_method="hist",
        scale_pos_weight=scale_pos_weight,
    )

    preds = np.zeros(len(rows), dtype=float)
    fold_summaries = []

    gkf = GroupKFold(n_splits=n_splits)

    for fold, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups), start=1):
        y_train = y[train_idx]
        y_test = y[test_idx]

        if len(set(y_train)) < 2:
            # This can happen with tiny grouped data.
            # Use a neutral probability for the held-out fold.
            preds[test_idx] = y_train.mean()
            fold_summaries.append({
                "fold": fold,
                "train_groups": sorted(set(groups[train_idx])),
                "test_groups": sorted(set(groups[test_idx])),
                "test_n": int(len(test_idx)),
                "test_pos": int(y_test.sum()),
                "test_neg": int((y_test == 0).sum()),
                "skipped_reason": "training fold had only one class",
            })
            continue

        model = XGBClassifier(**params)
        model.fit(X.iloc[train_idx], y_train)

        p = model.predict_proba(X.iloc[test_idx])[:, 1]
        preds[test_idx] = p
        yhat = (p >= 0.5).astype(int)

        fold_info = {
            "fold": fold,
            "train_groups": sorted(set(groups[train_idx])),
            "test_groups": sorted(set(groups[test_idx])),
            "test_n": int(len(test_idx)),
            "test_pos": int(y_test.sum()),
            "test_neg": int((y_test == 0).sum()),
            "accuracy_0.5": float(accuracy_score(y_test, yhat)),
        }

        if len(set(y_test)) == 2:
            fold_info["roc_auc"] = float(roc_auc_score(y_test, p))
            fold_info["average_precision"] = float(average_precision_score(y_test, p))
        else:
            fold_info["roc_auc"] = None
            fold_info["average_precision"] = None

        pr, rc, f1, _ = precision_recall_fscore_support(
            y_test, yhat, average="binary", zero_division=0
        )
        fold_info["precision_0.5"] = float(pr)
        fold_info["recall_0.5"] = float(rc)
        fold_info["f1_0.5"] = float(f1)

        fold_summaries.append(fold_info)

    final_model = XGBClassifier(**params)
    final_model.fit(X, y)

    Path(args.out_model_json).parent.mkdir(parents=True, exist_ok=True)
    final_model.save_model(args.out_model_json)

    yhat_all = (preds >= 0.5).astype(int)

    metrics = {
        "rows": int(len(rows)),
        "positives": int(pos),
        "negatives": int(neg),
        "positive_rate": float(pos / len(rows)),
        "source_windows": dict(source_counts),
        "physical_groups": dict(group_counts),
        "n_splits": int(n_splits),
        "features": FEATURES,
        "cv_accuracy_0.5": float(accuracy_score(y, yhat_all)),
        "folds": fold_summaries,
        "note": "Smoke model. chr20 30x/35x are collapsed into one physical group for CV.",
    }

    if len(set(y)) == 2:
        metrics["cv_roc_auc"] = float(roc_auc_score(y, preds))
        metrics["cv_average_precision"] = float(average_precision_score(y, preds))
    else:
        metrics["cv_roc_auc"] = None
        metrics["cv_average_precision"] = None

    Path(args.out_metrics_json).write_text(json.dumps(metrics, indent=2) + "\n")

    booster = final_model.get_booster()
    gain = booster.get_score(importance_type="gain")
    weight = booster.get_score(importance_type="weight")

    with open(args.out_feature_importance_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["feature", "gain", "weight"], delimiter="\t")
        w.writeheader()
        for f in FEATURES:
            w.writerow({
                "feature": f,
                "gain": gain.get(f, 0.0),
                "weight": weight.get(f, 0.0),
            })

    out_fields = list(rows[0].keys()) + [
        "physical_group",
        "xgb_oof_probability",
        "xgb_oof_pred_0.5",
    ]

    with open(args.out_predictions_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r, g, p, yh in zip(rows, groups, preds, yhat_all):
            rr = dict(r)
            rr["physical_group"] = g
            rr["xgb_oof_probability"] = f"{p:.6f}"
            rr["xgb_oof_pred_0.5"] = int(yh)
            w.writerow(rr)

    print(json.dumps(metrics, indent=2))

if __name__ == "__main__":
    main()
