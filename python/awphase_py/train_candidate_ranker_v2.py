#!/usr/bin/env python3
import argparse, csv, json, math
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def try_float(x):
    try: return float(x)
    except Exception: return math.nan

DROP_COLS = {
    'query_id','candidate_id','pos','candidate_type','truth_state','label',
    'panel_top_group','panel_top_group_local'
}

def query_top1_accuracy(rows, scores):
    best = {}
    for r, s in zip(rows, scores):
        q = r['query_id']
        if q not in best or s > best[q][0]:
            best[q] = (s, int(float(r.get('label',0) or 0)))
    return sum(v[1] for v in best.values()) / max(1, len(best))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--training-tsv', required=True)
    ap.add_argument('--out-model-json', required=True)
    ap.add_argument('--out-metrics-json', required=True)
    ap.add_argument('--out-feature-importance-tsv', required=True)
    args = ap.parse_args()

    rows = [r for r in read_tsv(args.training_tsv) if str(r.get('label','')).strip() != '']
    if not rows:
        raise SystemExit('No labeled rows in training TSV')

    cols = []
    for k in rows[0].keys():
        if k in DROP_COLS:
            continue
        vals = [try_float(r.get(k,'')) for r in rows[:50]]
        if any(not math.isnan(v) for v in vals):
            cols.append(k)

    X = []
    y = []
    qids = []
    weights = []
    for r in rows:
        y.append(int(float(r.get('label',0) or 0)))
        qids.append(r['query_id'])
        feat = []
        for c in cols:
            v = try_float(r.get(c,''))
            feat.append(0.0 if math.isnan(v) else v)
        X.append(feat)
        w = 1.0
        w += 0.5 * (0 if int(float(r.get('candidate_phase_state',0) or 0)) == 0 else 1)
        w += 0.5 * (1 if str(r.get('f_anchor_consistent','0')) == '1' else 0)
        w += 0.25 * (1 if str(r.get('f_candidate_matches_donor_local','0')) == '1' else 0)
        weights.append(w)

    import numpy as np, xgboost as xgb
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)

    unique_qids = []
    group_sizes = []
    cur = None
    count = 0
    for q in qids:
        if q != cur:
            if cur is not None:
                unique_qids.append(cur); group_sizes.append(count)
            cur = q; count = 1
        else:
            count += 1
    if cur is not None:
        unique_qids.append(cur); group_sizes.append(count)

    ranker = xgb.XGBRanker(
        objective='rank:pairwise',
        n_estimators=250,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.9,
        colsample_bytree=0.9,
        reg_lambda=1.0,
        min_child_weight=5,
        random_state=13,
    )
    ranker.fit(X, y, group=group_sizes, verbose=False)

    scores = ranker.predict(X)
    acc05 = sum((s >= 0.0) == bool(lbl) for s, lbl in zip(scores, y)) / max(1, len(y))
    qacc = query_top1_accuracy(rows, scores)

    ranker.save_model(args.out_model_json)
    meta = {'feature_names': cols, 'objective': 'rank:pairwise'}
    Path(args.out_model_json + '.meta.json').write_text(json.dumps(meta, indent=2))

    booster = ranker.get_booster()
    imp = booster.get_score(importance_type='gain')
    imp_rows = [{'feature': f, 'gain': float(imp.get(f, 0.0))} for f in cols]
    imp_rows.sort(key=lambda r: r['gain'], reverse=True)
    write_tsv(args.out_feature_importance_tsv, imp_rows, ['feature','gain'])

    metrics = {
        'rows': len(rows),
        'positive_rate': float(y.mean()),
        'accuracy_0.5': acc05,
        'query_top1_accuracy': qacc,
        'objective': 'rank:pairwise',
    }
    Path(args.out_metrics_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(metrics, open(args.out_metrics_json,'w'), indent=2)
    print(json.dumps(metrics, indent=2))

if __name__ == '__main__':
    main()
