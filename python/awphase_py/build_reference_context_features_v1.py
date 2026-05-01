#!/usr/bin/env python3
import argparse, csv, json, math
from pathlib import Path
import pysam


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)


def load_variants(path):
    obj = json.load(open(path))
    vals = obj.get('variants', obj) if isinstance(obj, dict) else obj
    rows = []
    for r in vals:
        try:
            rows.append((int(r['pos']), str(r.get('ref_allele', r.get('ref',''))), str(r.get('alt_allele', r.get('alt','')))))
        except Exception:
            pass
    rows.sort(key=lambda x: x[0])
    return rows


def shannon_entropy(seq):
    seq = [b for b in seq if b in 'ACGT']
    if not seq:
        return 0.0
    n = len(seq)
    e = 0.0
    for b in 'ACGT':
        c = seq.count(b)
        if c:
            p = c / n
            e -= p * math.log2(p)
    return e


def longest_hp(seq):
    best = cur = 0
    prev = None
    for b in seq:
        if b == prev:
            cur += 1
        else:
            cur = 1
        prev = b
        best = max(best, cur)
    return best


def side_hp(seq, center_idx, direction='left'):
    if not seq:
        return 0
    if center_idx < 0 or center_idx >= len(seq):
        return 0
    base = seq[center_idx]
    if base not in 'ACGT':
        return 0
    n = 0
    if direction == 'left':
        i = center_idx
        while i >= 0 and seq[i] == base:
            n += 1; i -= 1
    else:
        i = center_idx
        while i < len(seq) and seq[i] == base:
            n += 1; i += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--variant-json', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--flank', type=int, default=12)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    vars_ = load_variants(args.variant_json)
    fa = pysam.FastaFile(args.fasta)
    positions = [p for p,_,_ in vars_]
    rows = []
    for i, (pos, ref, alt) in enumerate(vars_):
        start = max(0, pos - 1 - args.flank)
        end = pos - 1 + max(len(ref),1) + args.flank
        seq = fa.fetch(args.chrom, start, end).upper()
        center = min(args.flank, len(seq)-1) if seq else 0
        gc = sum(1 for b in seq if b in 'GC') / max(1, sum(1 for b in seq if b in 'ACGT'))
        ent = shannon_entropy(seq)
        hp = longest_hp(seq) if seq else 0
        left_hp = side_hp(seq, center, 'left') if seq else 0
        right_hp = side_hp(seq, center, 'right') if seq else 0
        prev_dist = pos - positions[i-1] if i > 0 else 10**9
        next_dist = positions[i+1] - pos if i+1 < len(positions) else 10**9
        low_complex = 1 if ent < 1.5 or hp >= 6 else 0
        is_indel = 1 if len(ref) != len(alt) else 0
        rows.append({
            'pos': pos,
            'ctx_gc_frac': f'{gc:.6f}',
            'ctx_entropy': f'{ent:.6f}',
            'ctx_longest_hp': hp,
            'left_hp': left_hp,
            'right_hp': right_hp,
            'neighbor_dist_prev': prev_dist,
            'neighbor_dist_next': next_dist,
            'low_complexity_flag': low_complex,
            'is_indel': is_indel,
        })
    write_tsv(args.out_tsv, rows, list(rows[0].keys()) if rows else ['pos'])
    summ = {'rows': len(rows), 'chrom': args.chrom, 'flank': args.flank}
    json.dump(summ, open(args.out_summary_json, 'w'), indent=2)
    print(json.dumps(summ, indent=2))

if __name__ == '__main__':
    main()
