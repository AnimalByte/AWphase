#!/usr/bin/env python3
import argparse, csv, json, os
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--training-tsv', required=True)
    ap.add_argument('--out-model-json', required=True)
    ap.add_argument('--out-metrics-json', required=True)
    ap.add_argument('--out-feature-importance-tsv', required=True)
    args = ap.parse_args()

    import xgboost as xgb

    df = pd.read_csv(args.training_tsv, sep='\t')
    df = df[df['label'].notna() & (df['label'].astype(str) != '')].copy()
    df['label'] = df['label'].astype(int)
    # sort by query for ranking groups
    df.sort_values(['query_id', 'candidate_id'], inplace=True)
    exclude = {'query_id','candidate_id','candidate_type','label','truth_state'}
    feat_cols = [c for c in df.columns if c not in exclude]
    # keep only numeric features
    for c in feat_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0.0)
    X = df[feat_cols]
    y = df['label'].values
    group_sizes = df.groupby('query_id').size().tolist()

    model = xgb.XGBRanker(
        objective='rank:pairwise',
        eval_metric='map',
        n_estimators=300,
        learning_rate=0.05,
        max_depth=5,
        min_child_weight=4,
        subsample=0.85,
        colsample_bytree=0.85,
        reg_lambda=2.0,
        tree_method='hist',
        random_state=42,
    )
    model.fit(X, y, group=group_sizes, verbose=False)
    model.save_model(args.out_model_json)
    Path(args.out_model_json + '.features.json').write_text(json.dumps(feat_cols, indent=2))

    scores = model.predict(X)
    df['_score'] = scores
    acc05 = float(((scores >= 0.5).astype(int) == y).mean())
    top1_hits = []
    for _, g in df.groupby('query_id'):
        g2 = g.sort_values('_score', ascending=False)
        top1_hits.append(int(g2.iloc[0]['label'] == 1))
    top1 = float(np.mean(top1_hits)) if top1_hits else 0.0
    metrics = {
        'rows': int(len(df)),
        'positive_rate': float(df['label'].mean()) if len(df) else 0.0,
        'accuracy_0.5': acc05,
        'query_top1_accuracy': top1,
        'features': feat_cols,
        'objective': 'rank:pairwise',
    }
    json.dump(metrics, open(args.out_metrics_json, 'w'), indent=2)
    print(json.dumps({k: metrics[k] for k in ['rows','positive_rate','accuracy_0.5','query_top1_accuracy']}, indent=2))

    booster = model.get_booster()
    imp = booster.get_score(importance_type='gain')
    with open(args.out_feature_importance_tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['feature','gain'], delimiter='\t')
        w.writeheader()
        for f, g in sorted(imp.items(), key=lambda kv: kv[1], reverse=True):
            w.writerow({'feature': f, 'gain': g})

if __name__ == '__main__':
    main()
