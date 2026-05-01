#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def f(x, d=0.0):
    try: return float(x)
    except Exception: return d

def i(x, d=0):
    try: return int(float(x))
    except Exception: return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-scores-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--max-gap-bp', type=int, default=50000)
    ap.add_argument('--max-region-sites', type=int, default=8)
    ap.add_argument('--uncertain-margin', type=float, default=0.20)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    scaffold = {}
    for r in read_tsv(args.scaffold_calls_tsv):
        try:
            scaffold[int(r['pos'])] = int(float(r.get('scaffold_phase_state',0) or 0))
        except Exception:
            pass

    perq = {}
    for r in read_tsv(args.candidate_scores_tsv):
        perq.setdefault(r['query_id'], []).append(r)

    queries = []
    for qid, arr in perq.items():
        arr.sort(key=lambda r: f(r.get('model_score_v2', r.get('model_score', 0.0))), reverse=True)
        pos = i(arr[0]['pos'])
        best = arr[0]
        margin = f(best.get('query_margin_v2', best.get('query_margin', 0.0)))
        best_state = i(best.get('candidate_phase_state', 0))
        sc = scaffold.get(pos, 0)
        uncertain = (margin < args.uncertain_margin) or (sc == 0) or (best_state != 0 and sc != 0 and best_state != sc)
        queries.append({'query_id': qid, 'pos': pos, 'best_state': best_state, 'margin': margin, 'uncertain': 1 if uncertain else 0})

    queries.sort(key=lambda r: r['pos'])
    out = []
    region_id = 0
    cur_count = 0
    last_pos = None
    for q in queries:
        if q['uncertain'] != 1:
            continue
        pos = q['pos']
        if last_pos is None or (pos - last_pos) > args.max_gap_bp or cur_count >= args.max_region_sites:
            region_id += 1
            cur_count = 0
        out.append({
            'region_id': f'region_{region_id}',
            'query_id': q['query_id'],
            'pos': pos,
            'region_index': cur_count,
            'best_state': q['best_state'],
            'query_margin_v2': f'{q["margin"]:.6f}',
        })
        last_pos = pos
        cur_count += 1

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else ['region_id','query_id','pos','region_index','best_state','query_margin_v2'])
    summary = {'regions': len({r['region_id'] for r in out}), 'query_rows': len(out)}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
