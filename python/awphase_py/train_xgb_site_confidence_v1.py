#!/usr/bin/env python3
import argparse
import json
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--out-model", required=True)
    ap.add_argument("--out-metrics-json", required=True)
    ap.add_argument("--out-feature-importance-tsv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.training_tsv, sep="\t")

    # only phased calls for v1
    df = df[df["phase_state"] != 0].copy()
    df["y"] = (df["label"] == 1).astype(int)

    features = [
        "confidence",
        "alt_support",
        "ref_support",
        "read_bias",
        "donor_margin",
        "combined_donor_bias",
        "path_sign",
        "read_sign",
        "donor_sign",
        "path_margin_agrees",
        "donor_read_conflict",
        "weak_read_evidence",
        "forced_abstain",
        "neighbor_dist",
        "reads_seen",
        "ref_obs",
        "alt_obs",
        "other_base",
        "deletions",
        "refskips",
        "low_mapq",
        "low_baseq",
        "secondary_or_supplementary",
        "duplicates",
        "emitted_obs",
        "usable_site",
    ]

    X = df[features].fillna(0)
    y = df["y"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=42, stratify=y
    )

    model = xgb.XGBClassifier(
        n_estimators=400,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=42,
    )

    model.fit(X_train, y_train)

    prob = model.predict_proba(X_test)[:, 1]
    pred = (prob >= 0.5).astype(int)

    metrics = {
        "n_train": int(len(X_train)),
        "n_test": int(len(X_test)),
        "positive_rate_train": float(y_train.mean()),
        "positive_rate_test": float(y_test.mean()),
        "roc_auc": float(roc_auc_score(y_test, prob)),
        "average_precision": float(average_precision_score(y_test, prob)),
        "classification_report": classification_report(y_test, pred, output_dict=True),
        "features": features,
    }

    model.save_model(args.out_model)

    with open(args.out_metrics_json, "w") as fh:
        json.dump(metrics, fh, indent=2)

    imp = pd.DataFrame({
        "feature": features,
        "importance": model.feature_importances_,
    }).sort_values("importance", ascending=False)
    imp.to_csv(args.out_feature_importance_tsv, sep="\t", index=False)

    print(json.dumps({
        "roc_auc": metrics["roc_auc"],
        "average_precision": metrics["average_precision"],
        "n_train": metrics["n_train"],
        "n_test": metrics["n_test"],
    }, indent=2))

if __name__ == "__main__":
    main()
