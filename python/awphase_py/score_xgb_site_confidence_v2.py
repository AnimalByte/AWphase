#!/usr/bin/env python3
import argparse
import json
import pandas as pd
import xgboost as xgb

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inference-tsv", required=True)
    ap.add_argument("--model-json", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.inference_tsv, sep="\t")

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

    model = xgb.XGBClassifier()
    model.load_model(args.model_json)

    booster = model.get_booster()
    feat_names = booster.feature_names
    for f in feat_names:
        if f not in X.columns:
            X[f] = 0
    X = X[feat_names]

    prob = model.predict_proba(X)[:, 1]

    out = df.copy()
    out["xgb_safe_prob"] = prob
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(out[["pos", "phase_state", "xgb_safe_prob"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
