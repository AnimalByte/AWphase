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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inference-tsv', required=True)
    ap.add_argument('--model-json', required=True)
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()

    rows = read_tsv(args.inference_tsv)
    meta = json.load(open(args.model_json + '.meta.json'))
    cols = meta['feature_names']

    import numpy as np, xgboost as xgb
    X = []
    for r in rows:
        feat = []
        for c in cols:
            v = try_float(r.get(c,''))
            feat.append(0.0 if math.isnan(v) else v)
        X.append(feat)
    X = np.asarray(X, dtype=float)

    model = xgb.XGBRanker()
    model.load_model(args.model_json)
    scores = model.predict(X)

    byq = {}
    for r, s in zip(rows, scores):
        q = r['query_id']
        byq.setdefault(q, []).append((s, r))
    out = []
    for q, arr in byq.items():
        arr.sort(key=lambda x: x[0], reverse=True)
        best = arr[0][0]
        second = arr[1][0] if len(arr) > 1 else arr[0][0]
        for rank, (s, r) in enumerate(arr, start=1):
            row = dict(r)
            row['model_score_v2'] = f'{float(s):.6f}'
            row['query_rank_v2'] = rank
            row['query_best_score_v2'] = f'{float(best):.6f}'
            row['query_margin_v2'] = f'{float(best-second):.6f}'
            out.append(row)

    fields = list(out[0].keys()) if out else []
    write_tsv(args.out_tsv, out, fields)
    print(f'Wrote {len(out)} rows to {args.out_tsv}')

if __name__ == '__main__':
    main()
