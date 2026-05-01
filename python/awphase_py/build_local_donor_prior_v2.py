#!/usr/bin/env python3
import argparse, csv, json, math
from collections import Counter
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def _looks_like_int(x):
    try: int(str(x).strip()); return True
    except Exception: return False

def load_target_positions(path):
    positions = []
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh if ln.strip()]
    if not lines:
        raise SystemExit(f'Empty target sites file: {path}')
    first = lines[0].split('\t')
    if len(first) >= 2 and not _looks_like_int(first[0]) and _looks_like_int(first[1]):
        for ln in lines:
            p = ln.split('\t')
            if len(p) >= 2 and _looks_like_int(p[1]):
                positions.append(int(p[1]))
    elif len(first) >= 1 and _looks_like_int(first[0]):
        for ln in lines:
            p = ln.split('\t')
            if _looks_like_int(p[0]):
                positions.append(int(p[0]))
    else:
        rows = read_tsv(path)
        if not rows: raise SystemExit(f'No rows in {path}')
        pos_col = None
        for c in rows[0].keys():
            lc = str(c).lower()
            if lc in {'pos','position','site_pos'} or 'pos' in lc:
                pos_col = c; break
        if pos_col is None:
            raise SystemExit(f'Could not find pos column in {path}')
        for r in rows:
            try: positions.append(int(r[pos_col]))
            except Exception: pass
    return positions

def load_scaffold(path):
    states = {}
    for r in read_tsv(path):
        try:
            pos = int(r['pos']); state = int(float(r.get('scaffold_phase_state',0) or 0))
            conf = float(r.get('scaffold_confidence',0.0) or 0.0)
        except Exception:
            continue
        if state != 0:
            states[pos] = (state, conf)
    return states

def load_panel_candidates(path):
    obj = json.load(open(path))
    if not isinstance(obj, list):
        raise SystemExit('panel candidates expected list')
    out = []
    for rec in obj:
        alleles = rec.get('alleles')
        if not isinstance(alleles, list):
            continue
        out.append({
            'donor_id': str(rec.get('donor_id','')),
            'hap_id': rec.get('hap_id',''),
            'group': str(rec.get('group','')),
            'alleles': alleles,
        })
    if not out:
        raise SystemExit(f'No usable donor records in {path}')
    return out

def entropy_binary(p):
    if p <= 0.0 or p >= 1.0: return 0.0
    return -(p*math.log2(p) + (1-p)*math.log2(1-p))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-queries-tsv', required=True)
    ap.add_argument('--panel-candidates-json', required=True)
    ap.add_argument('--target-sites-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--window-sites', type=int, default=8)
    ap.add_argument('--top-donors', type=int, default=64)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    positions = load_target_positions(args.target_sites_tsv)
    pos_to_idx = {p:i for i,p in enumerate(positions)}
    scaffold = load_scaffold(args.scaffold_calls_tsv)
    donors = load_panel_candidates(args.panel_candidates_json)
    queries = read_tsv(args.candidate_queries_tsv)

    out = []
    total_support = 0.0
    pos_hits = 0
    groups = Counter()
    cache = {}

    for q in queries:
        pos = int(q['pos'])
        state = int(q['candidate_phase_state'])
        idx = pos_to_idx.get(pos)
        top_group = ''
        donor_bias = 0.0
        donor_entropy = 1.0
        panel_support_sum = 0.0
        panel_support_frac = 0.0
        panel_oppose_sum = 0.0
        panel_total_informative = 0.0
        candidate_matches_donor = 0
        donor_window_anchor_support = 0.0
        donor_window_anchor_disagree = 0.0
        donor_window_entropy = 1.0
        top_donor_mean_score = 0.0

        if idx is not None and state != 0:
            if idx not in cache:
                lo = max(0, idx - args.window_sites)
                hi = min(len(positions)-1, idx + args.window_sites)
                anchor_idxs = []
                for j in range(lo, hi+1):
                    p = positions[j]
                    if p in scaffold:
                        anchor_idxs.append((j, scaffold[p][0], scaffold[p][1]))
                donor_scores = []
                for d in donors:
                    arr = d['alleles']
                    score = 0.0; agree = 0.0; disagree = 0.0
                    for j,s,c in anchor_idxs:
                        if j >= len(arr): continue
                        a = arr[j]
                        if a == 0: continue
                        w = max(0.1, c)
                        if a == s:
                            score += w; agree += w
                        elif a == -s:
                            score -= w; disagree += w
                    donor_scores.append((score, agree, disagree, d))
                donor_scores.sort(key=lambda x: (x[0], x[1]), reverse=True)
                cache[idx] = donor_scores[:args.top_donors]

            top = cache[idx]
            support = oppose = 0.0
            top_score_sum = 0.0
            grp_counter = Counter()
            agree_total = disagree_total = 0.0
            for score, agree, disagree, d in top:
                arr = d['alleles']
                if idx >= len(arr): continue
                a = arr[idx]
                top_score_sum += score
                agree_total += agree
                disagree_total += disagree
                if a == 0: 
                    continue
                grp_counter[d['group']] += 1
                if a == state:
                    support += 1.0
                elif a == -state:
                    oppose += 1.0

            informative = support + oppose
            panel_support_sum = support
            panel_oppose_sum = oppose
            panel_total_informative = informative
            top_donor_mean_score = top_score_sum / max(1, len(top))
            donor_window_anchor_support = agree_total
            donor_window_anchor_disagree = disagree_total
            if informative > 0:
                pos_hits += 1
                panel_support_frac = support / informative
                donor_bias = (support - oppose) / informative
                donor_entropy = entropy_binary(panel_support_frac)
                candidate_matches_donor = 1 if support > oppose else 0
                denom = agree_total + disagree_total
                donor_window_entropy = entropy_binary(agree_total / denom) if denom > 0 else 1.0
                if grp_counter:
                    top_group = grp_counter.most_common(1)[0][0]
                    groups[top_group] += 1
            total_support += panel_support_sum

        out.append({
            'query_id': q['query_id'],
            'candidate_id': q['candidate_id'],
            'pos': pos,
            'candidate_phase_state': state,
            'donor_bias_local': f'{donor_bias:.6f}',
            'donor_entropy_local': f'{donor_entropy:.6f}',
            'candidate_matches_donor_local': candidate_matches_donor,
            'panel_support_sum_local': f'{panel_support_sum:.6f}',
            'panel_support_frac_local': f'{panel_support_frac:.6f}',
            'panel_oppose_sum_local': f'{panel_oppose_sum:.6f}',
            'panel_total_informative_local': f'{panel_total_informative:.6f}',
            'panel_top_group_local': top_group,
            'top_donor_mean_score_local': f'{top_donor_mean_score:.6f}',
            'donor_window_anchor_support': f'{donor_window_anchor_support:.6f}',
            'donor_window_anchor_disagree': f'{donor_window_anchor_disagree:.6f}',
            'donor_window_entropy': f'{donor_window_entropy:.6f}',
        })

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else [])
    summary = {
        'rows': len(out),
        'panel_top_group': groups.most_common(1)[0][0] if groups else '',
        'panel_support_sum': total_support,
        'panel_pos_specific_hits': pos_hits,
        'panel_positions': len(positions),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
