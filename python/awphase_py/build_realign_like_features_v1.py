#!/usr/bin/env python3
import argparse, csv, json
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def load_haplotag(path):
    tags = {}
    for r in read_tsv(path):
        rid = str(r.get('read_id','')).strip()
        hp = str(r.get('hp', r.get('HP', r.get('haplotype','UNK')))).upper()
        if hp in {'1','HP1'}: hp = 'HP1'
        elif hp in {'2','HP2'}: hp = 'HP2'
        else: hp = 'UNK'
        tags[rid] = hp
    return tags

def load_obs_by_pos(path, hap_tags):
    by = defaultdict(list)
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            rid = str(r.get('read_id','')).strip()
            hp = hap_tags.get(rid, 'UNK')
            try:
                pos = int(r.get('site_pos', r.get('pos')))
                allele = int(float(r.get('allele', 0) or 0))
                delta = float(r.get('allele_score_delta', 0.0) or 0.0)
                conf = float(r.get('allele_confidence', 0.0) or 0.0)
                ambig = int(float(r.get('is_ambiguous', 0) or 0))
            except Exception:
                continue
            if ambig == 1 or allele not in (-1, 1):
                continue
            by[pos].append((hp, allele, delta, conf))
    return by

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-queries-tsv', required=True)
    ap.add_argument('--obs-tsv', required=True)
    ap.add_argument('--haplotag-tsv', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    hap = load_haplotag(args.haplotag_tsv)
    by = load_obs_by_pos(args.obs_tsv, hap)
    out = []

    for q in read_tsv(args.candidate_queries_tsv):
        pos = int(q['pos']); state = int(q['candidate_phase_state'])
        support = oppose = 0.0
        sup_delta = opp_delta = 0.0
        tagged_reads = 0
        for hp, allele, delta, conf in by.get(pos, []):
            w = max(0.01, abs(delta) * max(0.01, conf))
            if hp in {'HP1','HP2'}:
                tagged_reads += 1
                if state == 1:
                    ok = (hp == 'HP1' and allele == 1) or (hp == 'HP2' and allele == -1)
                elif state == -1:
                    ok = (hp == 'HP1' and allele == -1) or (hp == 'HP2' and allele == 1)
                else:
                    ok = False
                if ok:
                    support += w; sup_delta += abs(delta)
                else:
                    oppose += w; opp_delta += abs(delta)
        total = support + oppose
        margin = (support - oppose) / total if total > 0 else 0.0
        out.append({
            'query_id': q['query_id'],
            'candidate_id': q['candidate_id'],
            'pos': pos,
            'candidate_phase_state': state,
            'realign_support_w': f'{support:.6f}',
            'realign_oppose_w': f'{oppose:.6f}',
            'realign_margin': f'{margin:.6f}',
            'realign_support_delta_sum': f'{sup_delta:.6f}',
            'realign_oppose_delta_sum': f'{opp_delta:.6f}',
            'realign_tagged_obs': tagged_reads,
        })

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else [])
    summary = {'rows': len(out)}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
