#!/usr/bin/env python3
import argparse, csv, itertools, json
from collections import defaultdict
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

def load_obs(obs_tsv, region_pos_set):
    by_read = defaultdict(list)
    with open(obs_tsv) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            try:
                pos = int(r.get('site_pos', r.get('pos')))
                allele = int(float(r.get('allele',0) or 0))
            except Exception:
                continue
            if pos not in region_pos_set or allele not in (-1, 1):
                continue
            w = f(r.get('weight_v2', 0.0), 0.0)
            if w <= 0:
                w = max(0.01, f(r.get('allele_confidence',0.0),0.0) + abs(f(r.get('allele_score_delta',0.0),0.0))*0.1)
            by_read[str(r.get('read_id',''))].append((pos, allele, w, 'read'))
    return by_read

def load_scaffold_reads(tsv, region_pos_set):
    by_read = defaultdict(list)
    for r in read_tsv(tsv):
        pos = i(r.get('site_pos', r.get('pos')))
        if pos not in region_pos_set: continue
        allele = i(r.get('allele', 0))
        if allele not in (-1, 1): continue
        w = f(r.get('weight', 0.0), 0.0)
        by_read[str(r.get('read_id',''))].append((pos, allele, w, 'scaffold'))
    return by_read

def pairwise_score(read_items, assign):
    items = [(p,a,w,src) for p,a,w,src in read_items if assign.get(p,0) != 0]
    if len(items) < 2:
        return 0.0
    score = 0.0
    for j in range(len(items)):
        p1,a1,w1,src1 = items[j]
        for k in range(j+1, len(items)):
            p2,a2,w2,src2 = items[k]
            rel_obs = a1 * a2
            rel_assign = assign[p1] * assign[p2]
            w = min(w1, w2)
            if src1 == 'scaffold' or src2 == 'scaffold':
                w *= 1.5
            score += w if rel_obs == rel_assign else -w
    return score

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--region-queries-tsv', required=True)
    ap.add_argument('--candidate-scores-tsv', required=True)
    ap.add_argument('--obs-tsv', required=True)
    ap.add_argument('--scaffold-as-reads-tsv', required=True)
    ap.add_argument('--max-exact-sites', type=int, default=8)
    ap.add_argument('--beam-size', type=int, default=64)
    ap.add_argument('--weight-score', type=float, default=1.0)
    ap.add_argument('--weight-pairwise', type=float, default=0.35)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    region_rows = read_tsv(args.region_queries_tsv)
    cand_rows = read_tsv(args.candidate_scores_tsv)

    per_region = defaultdict(list)
    query_to_region = {}
    for r in region_rows:
        rid = r['region_id']; q = r['query_id']
        per_region[rid].append(r)
        query_to_region[q] = rid

    cand_by_query = defaultdict(list)
    for r in cand_rows:
        q = r['query_id']
        if q in query_to_region:
            cand_by_query[q].append(r)

    out = []
    chosen_regions = 0
    resolved_positions = 0

    for rid, qrows in per_region.items():
        qrows = sorted(qrows, key=lambda r: i(r['pos']))
        positions = [i(r['pos']) for r in qrows]
        region_pos_set = set(positions)
        obs_reads = load_obs(args.obs_tsv, region_pos_set)
        scaf_reads = load_scaffold_reads(args.scaffold_as_reads_tsv, region_pos_set)
        read_sets = list(obs_reads.values()) + list(scaf_reads.values())

        sites = []
        for qr in qrows:
            q = qr['query_id']
            arr = sorted(cand_by_query[q], key=lambda r: f(r.get('model_score_v2', r.get('model_score', 0.0))), reverse=True)
            keep = []
            abstain = None
            for r in arr:
                st = i(r.get('candidate_phase_state',0))
                if st == 0:
                    abstain = r
                else:
                    keep.append(r)
            keep = keep[:2]
            if abstain is not None:
                keep.append(abstain)
            if not keep:
                keep = arr[:1]
            sites.append((q, i(qr['pos']), keep))

        def score_assignment(assignment_rows):
            assign = {pos: i(row.get('candidate_phase_state',0)) for _,pos,row in assignment_rows}
            local = 0.0
            for _, pos, row in assignment_rows:
                s = f(row.get('model_score_v2', row.get('model_score', 0.0)))
                s += 0.8 * f(row.get('f_anchor_consistent', 0.0))
                s -= 0.6 * f(row.get('f_anchor_conflict', 0.0))
                s += 0.4 * f(row.get('f_candidate_matches_donor_local', 0.0))
                s += 0.4 * f(row.get('f_realign_margin', 0.0))
                s += 0.2 * f(row.get('f_matches_left_anchor', 0.0))
                s += 0.2 * f(row.get('f_matches_right_anchor', 0.0))
                if i(row.get('candidate_phase_state',0)) == 0:
                    s -= 0.05
                local += s
            pair = sum(pairwise_score(rs, assign) for rs in read_sets)
            return args.weight_score * local + args.weight_pairwise * pair

        n = len(sites)
        total_states = 1
        for _,_,opts in sites: total_states *= len(opts)

        best = None
        second = None
        if n <= args.max_exact_sites and total_states <= 50000:
            for combo in itertools.product(*[opts for _,_,opts in sites]):
                ass = [(sites[idx][0], sites[idx][1], combo[idx]) for idx in range(n)]
                sc = score_assignment(ass)
                if best is None or sc > best[0]:
                    second = best
                    best = (sc, ass)
                elif second is None or sc > second[0]:
                    second = (sc, ass)
        else:
            beam = [(0.0, [])]
            for q,pos,opts in sites:
                next_beam = []
                for _, base_ass in beam:
                    for row in opts:
                        ass = base_ass + [(q,pos,row)]
                        sc = score_assignment(ass)
                        next_beam.append((sc, ass))
                next_beam.sort(key=lambda x: x[0], reverse=True)
                beam = next_beam[:args.beam_size]
            best = beam[0] if beam else (0.0, [])
            second = beam[1] if len(beam) > 1 else best

        best_score = best[0]
        second_score = second[0] if second else best_score
        margin = best_score - second_score
        chosen_regions += 1

        for q, pos, row in best[1]:
            st = i(row.get('candidate_phase_state',0))
            out.append({
                'region_id': rid,
                'query_id': q,
                'pos': pos,
                'resolved_phase_state': st,
                'resolved_score': f'{best_score:.6f}',
                'resolved_margin': f'{margin:.6f}',
                'query_best_score_v2': row.get('query_best_score_v2', row.get('query_best_score','0.0')),
                'query_margin_v2': row.get('query_margin_v2', row.get('query_margin','0.0')),
                'resolved_candidate_id': row.get('candidate_id',''),
                'resolved_anchor_consistent': row.get('f_anchor_consistent','0'),
                'resolved_donor_match': row.get('f_candidate_matches_donor_local','0'),
            })
            resolved_positions += 1

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else ['region_id','query_id','pos','resolved_phase_state'])
    summary = {'regions': chosen_regions, 'resolved_positions': resolved_positions}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
