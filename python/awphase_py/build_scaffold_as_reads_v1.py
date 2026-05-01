#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--min-confidence', type=float, default=0.65)
    ap.add_argument('--weight-scale', type=float, default=8.0)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    rows_in = read_tsv(args.scaffold_calls_tsv)
    out = []
    kept_blocks = set()
    for r in rows_in:
        bid = str(r.get('scaffold_block_id','')).strip()
        try:
            pos = int(r.get('pos'))
            state = int(float(r.get('scaffold_phase_state', 0) or 0))
            conf = float(r.get('scaffold_confidence', 0.0) or 0.0)
        except Exception:
            continue
        if not bid or state == 0 or conf < args.min_confidence:
            continue
        kept_blocks.add(bid)
        out.append({
            'read_id': f'SCAFFOLD::{bid}',
            'site_pos': pos,
            'allele': state,
            'weight': f"{max(0.1, conf * args.weight_scale):.6f}",
            'source': 'scaffold',
            'block_id': bid,
            'confidence': f"{conf:.6f}",
        })

    write_tsv(args.out_tsv, out, ['read_id','site_pos','allele','weight','source','block_id','confidence'])
    summary = {'rows': len(out), 'blocks': len(kept_blocks), 'min_confidence': args.min_confidence}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
