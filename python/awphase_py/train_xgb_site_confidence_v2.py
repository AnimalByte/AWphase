#!/usr/bin/env python3
import argparse
import json
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report
from sklearn.model_selection import GroupShuffleSplit

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--out-model", required=True)
    ap.add_argument("--out-metrics-json", required=True)
    ap.add_argument("--out-feature-importance-tsv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.training_tsv, sep="\t")

    # phased calls only for v2 confidence model
    df = df[df["phase_state"] != 0].copy()
    df["y"] = (df["label"] == 1).astype(int)

    numeric_features = [
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
        "ambiguous_obs",
        "alt_like_context_obs",
        "ref_like_context_obs",
        "mean_allele_confidence",
        "mean_signed_score_delta",
        "mean_abs_score_delta",
        "mean_local_ref_score",
        "mean_local_alt_score",
        "usable_site",
        "window_index",
    ]

    cat_features = ["decision_source"]

    X_num = df[numeric_features].fillna(0)
    X_cat = pd.get_dummies(df[cat_features].fillna("NA"), columns=cat_features)
    X = pd.concat([X_num, X_cat], axis=1)
    y = df["y"].astype(int)
    groups = df["group_id"]

    gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=42)
    train_idx, test_idx = next(gss.split(X, y, groups=groups))

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    model = xgb.XGBClassifier(
        n_estimators=600,
        max_depth=5,
        learning_rate=0.03,
        subsample=0.85,
        colsample_bytree=0.85,
        reg_lambda=1.0,
        min_child_weight=5,
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
        "n_groups_train": int(pd.Series(groups.iloc[train_idx]).nunique()),
        "n_groups_test": int(pd.Series(groups.iloc[test_idx]).nunique()),
        "numeric_features": numeric_features,
        "categorical_features": cat_features,
    }

    model.save_model(args.out_model)

    with open(args.out_metrics_json, "w") as fh:
        json.dump(metrics, fh, indent=2)

    imp = pd.DataFrame({
        "feature": X.columns,
        "importance": model.feature_importances_,
    }).sort_values("importance", ascending=False)
    imp.to_csv(args.out_feature_importance_tsv, sep="\t", index=False)

    print(json.dumps({
        "roc_auc": metrics["roc_auc"],
        "average_precision": metrics["average_precision"],
        "n_train": metrics["n_train"],
        "n_test": metrics["n_test"],
        "n_groups_train": metrics["n_groups_train"],
        "n_groups_test": metrics["n_groups_test"],
    }, indent=2))

if __name__ == "__main__":
    main()
