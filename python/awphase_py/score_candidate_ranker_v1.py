#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import pandas as pd
import xgboost as xgb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inference-tsv', required=True)
    ap.add_argument('--model-json', required=True)
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.inference_tsv, sep='\t')
    booster = xgb.Booster()
    booster.load_model(args.model_json)
    feat_path = args.model_json + '.features.json'
    if Path(feat_path).exists():
        feat_cols = json.load(open(feat_path))
    else:
        feat_cols = booster.feature_names
    for c in feat_cols:
        if c not in df.columns:
            df[c] = 0.0
    X = df[feat_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    dmat = xgb.DMatrix(X, feature_names=feat_cols)
    df['pred_score'] = booster.predict(dmat)

    # per-query ranking stats
    best_score = {}
    second_score = {}
    best_candidate = {}
    for qid, g in df.groupby('query_id'):
        vals = sorted(g['pred_score'].tolist(), reverse=True)
        best_score[qid] = vals[0] if vals else 0.0
        second_score[qid] = vals[1] if len(vals) > 1 else (vals[0] if vals else 0.0)
        best_candidate[qid] = g.sort_values('pred_score', ascending=False).iloc[0]['candidate_id'] if len(g) else ''
    df['query_best_score'] = df['query_id'].map(best_score)
    df['query_second_score'] = df['query_id'].map(second_score)
    df['query_margin'] = df['query_best_score'] - df['query_second_score']
    df['query_is_top'] = (df['candidate_id'] == df['query_id'].map(best_candidate)).astype(int)
    df.to_csv(args.out_tsv, sep='\t', index=False)
    print(f"Wrote {len(df)} rows to {args.out_tsv}")

if __name__ == '__main__':
    main()
