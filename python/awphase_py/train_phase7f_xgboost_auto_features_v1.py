#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, precision_recall_fscore_support
from sklearn.model_selection import GroupKFold
from xgboost import XGBClassifier

EXCLUDE = {
    "label",
    "source_window",
    "pos",
    "block_id",
    "pred_state",
    "truth_state",
    "block_orientation",
    "oriented_pred_state",
}

def physical_group(s):
    s = str(s)
    s = s.replace("_illumina30x_phase7a", "")
    s = s.replace("_illumina35x_phase7a", "")
    s = s.replace("_phase7a_general_runner", "")
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--out-model-json", required=True)
    ap.add_argument("--out-metrics-json", required=True)
    ap.add_argument("--out-feature-importance-tsv", required=True)
    ap.add_argument("--out-predictions-tsv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.training_tsv, sep="\t")

    if "label" not in df.columns or "source_window" not in df.columns:
        raise SystemExit("Need label and source_window columns")

    y = df["label"].astype(int).values
    groups = df["source_window"].map(physical_group).values

    candidate_features = []
    for col in df.columns:
        if col in EXCLUDE:
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        if numeric.notna().sum() == 0:
            continue
        if numeric.nunique(dropna=True) <= 1:
            continue
        candidate_features.append(col)

    X = pd.DataFrame({
        col: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        for col in candidate_features
    })

    unique_groups = sorted(set(groups))
    n_splits = min(5, len(unique_groups))

    pos = int((y == 1).sum())
    neg = int((y == 0).sum())

    params = dict(
        n_estimators=250,
        max_depth=3,
        learning_rate=0.035,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=2,
        reg_lambda=2.0,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=42,
        n_jobs=4,
        tree_method="hist",
        scale_pos_weight=(neg / pos if pos else 1.0),
    )

    preds = np.zeros(len(df), dtype=float)
    folds = []

    gkf = GroupKFold(n_splits=n_splits)

    for fold, (tr, te) in enumerate(gkf.split(X, y, groups), start=1):
        model = XGBClassifier(**params)
        model.fit(X.iloc[tr], y[tr])
        p = model.predict_proba(X.iloc[te])[:, 1]
        preds[te] = p

        yhat = (p >= 0.5).astype(int)
        info = {
            "fold": fold,
            "train_groups": sorted(set(groups[tr])),
            "test_groups": sorted(set(groups[te])),
            "test_n": int(len(te)),
            "test_pos": int(y[te].sum()),
            "test_neg": int((y[te] == 0).sum()),
            "accuracy_0.5": float(accuracy_score(y[te], yhat)),
        }

        if len(set(y[te])) == 2:
            info["roc_auc"] = float(roc_auc_score(y[te], p))
            info["average_precision"] = float(average_precision_score(y[te], p))
        else:
            info["roc_auc"] = None
            info["average_precision"] = None

        pr, rc, f1, _ = precision_recall_fscore_support(y[te], yhat, average="binary", zero_division=0)
        info["precision_0.5"] = float(pr)
        info["recall_0.5"] = float(rc)
        info["f1_0.5"] = float(f1)
        folds.append(info)

    final = XGBClassifier(**params)
    final.fit(X, y)

    Path(args.out_model_json).parent.mkdir(parents=True, exist_ok=True)
    final.save_model(args.out_model_json)

    yhat_all = (preds >= 0.5).astype(int)

    metrics = {
        "rows": int(len(df)),
        "positives": pos,
        "negatives": neg,
        "positive_rate": float(pos / len(df)),
        "source_windows": dict(Counter(df["source_window"])),
        "physical_groups": dict(Counter(groups)),
        "n_splits": n_splits,
        "n_features": len(candidate_features),
        "features": candidate_features,
        "cv_accuracy_0.5": float(accuracy_score(y, yhat_all)),
        "cv_roc_auc": float(roc_auc_score(y, preds)),
        "cv_average_precision": float(average_precision_score(y, preds)),
        "folds": folds,
    }

    Path(args.out_metrics_json).write_text(json.dumps(metrics, indent=2) + "\n")

    booster = final.get_booster()
    gain = booster.get_score(importance_type="gain")
    weight = booster.get_score(importance_type="weight")

    with open(args.out_feature_importance_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["feature", "gain", "weight"], delimiter="\t")
        w.writeheader()
        for f in candidate_features:
            w.writerow({"feature": f, "gain": gain.get(f, 0.0), "weight": weight.get(f, 0.0)})

    out_df = df.copy()
    out_df["physical_group"] = groups
    out_df["xgb_oof_probability"] = preds
    out_df["xgb_oof_pred_0.5"] = yhat_all
    out_df.to_csv(args.out_predictions_tsv, sep="\t", index=False)

    print(json.dumps(metrics, indent=2))

if __name__ == "__main__":
    main()
