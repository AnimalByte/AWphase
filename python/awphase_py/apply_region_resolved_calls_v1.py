#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    out_fields = list(fields) if fields else []
    seen = set(out_fields)
    for r in rows:
        for k in r.keys():
            if k not in seen:
                out_fields.append(k); seen.add(k)
    with open(path,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def f(x,d=0.0):
    try: return float(x)
    except Exception: return d
def i(x,d=0):
    try: return int(float(x))
    except Exception: return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--local-calls-tsv', required=True)
    ap.add_argument('--region-resolutions-tsv', required=True)
    ap.add_argument('--apply-margin', type=float, default=0.75)
    ap.add_argument('--min-query-margin', type=float, default=0.15)
    ap.add_argument('--require-anchor-or-donor', action='store_true')
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()

    res = {}
    for r in read_tsv(args.region_resolutions_tsv):
        pos = i(r['pos'])
        res[pos] = r

    rows = read_tsv(args.local_calls_tsv)
    out = []
    modified = 0
    for r in rows:
        row = dict(r)
        pos = i(r.get('pos',0))
        rr = res.get(pos)
        applied = 0
        if rr is not None:
            resolved = i(rr.get('resolved_phase_state',0))
            margin = f(rr.get('resolved_margin',0.0))
            qmargin = f(rr.get('query_margin_v2',0.0))
            anchor = i(rr.get('resolved_anchor_consistent',0))
            donor = i(rr.get('resolved_donor_match',0))
            ok = margin >= args.apply_margin and qmargin >= args.min_query_margin
            if args.require_anchor_or_donor:
                ok = ok and (anchor == 1 or donor == 1)
            if ok and resolved != 0:
                cur = i(r.get('local_phase_state', r.get('phase_state',0)))
                if cur != resolved:
                    if 'local_phase_state' in row:
                        row['local_phase_state'] = resolved
                    else:
                        row['phase_state'] = resolved
                    applied = 1
                    modified += 1
            row['phase4_resolved_score'] = rr.get('resolved_score','')
            row['phase4_resolved_margin'] = rr.get('resolved_margin','')
            row['phase4_query_margin'] = rr.get('query_margin_v2','')
            row['phase4_anchor_consistent'] = rr.get('resolved_anchor_consistent','')
            row['phase4_donor_match'] = rr.get('resolved_donor_match','')
        row['phase4_applied'] = applied
        out.append(row)

    write_tsv(args.out_tsv, out, list(rows[0].keys()) if rows else None)
    print({'rows': len(out), 'modified_rows': modified, 'apply_margin': args.apply_margin})

if __name__ == '__main__':
    main()
