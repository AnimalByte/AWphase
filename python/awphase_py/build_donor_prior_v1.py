#!/usr/bin/env python3
import argparse, csv, json, math
from collections import defaultdict, Counter
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def _looks_like_int(x):
    try:
        int(str(x).strip())
        return True
    except Exception:
        return False


def load_target_positions(path):
    positions = []
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh if ln.strip()]
    if not lines:
        raise SystemExit(f"Target-sites TSV is empty: {path}")
    first = lines[0].split('\t')
    if len(first) >= 2 and not _looks_like_int(first[0]) and _looks_like_int(first[1]):
        for ln in lines:
            parts = ln.split('\t')
            if len(parts) >= 2 and _looks_like_int(parts[1]):
                positions.append(int(parts[1]))
    elif len(first) >= 1 and _looks_like_int(first[0]):
        for ln in lines:
            parts = ln.split('\t')
            if len(parts) >= 1 and _looks_like_int(parts[0]):
                positions.append(int(parts[0]))
    else:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not rows:
            raise SystemExit(f"Target-sites TSV has no data rows: {path}")
        pos_col = None
        for c in rows[0].keys():
            lc = str(c).lower()
            if lc in {'pos', 'position', 'site_pos'} or 'pos' in lc:
                pos_col = c
                break
        if pos_col is None:
            raise SystemExit(f"Could not find position column in target-sites TSV. Columns={list(rows[0].keys())}")
        for r in rows:
            try:
                positions.append(int(r[pos_col]))
            except Exception:
                pass
    if not positions:
        raise SystemExit(f"No usable positions loaded from target-sites TSV: {path}")
    return positions


def load_panel_candidates(path):
    obj = json.load(open(path))
    if not isinstance(obj, list):
        raise SystemExit('panel_candidates JSON is expected to be a list of donor haplotypes')
    out = []
    for rec in obj:
        if not isinstance(rec, dict):
            continue
        alleles = rec.get('alleles')
        if not isinstance(alleles, list):
            continue
        out.append({
            'donor_id': str(rec.get('donor_id', '')),
            'hap_id': rec.get('hap_id', ''),
            'group': str(rec.get('group', '')),
            'alleles': alleles,
        })
    if not out:
        raise SystemExit(f'No usable donor records found in {path}')
    return out


def load_scaffold(scaffold_calls_tsv):
    rows = read_tsv(scaffold_calls_tsv)
    out = {}
    for r in rows:
        try:
            pos = int(r['pos'])
            state = int(float(r.get('scaffold_phase_state', 0) or 0))
            conf = float(r.get('scaffold_confidence', 0.0) or 0.0)
        except Exception:
            continue
        if state != 0:
            out[pos] = (state, conf)
    return out


def softmax(xs):
    if not xs:
        return []
    m = max(xs)
    ex = [math.exp(x - m) for x in xs]
    s = sum(ex)
    return [v / s for v in ex] if s > 0 else [1.0 / len(xs)] * len(xs)


def entropy(ps):
    e = 0.0
    for p in ps:
        if p > 0:
            e -= p * math.log2(p)
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-queries-tsv', required=True)
    ap.add_argument('--panel-candidates-json', required=True)
    ap.add_argument('--target-sites-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--window-sites', type=int, default=12)
    ap.add_argument('--min-anchor-count', type=int, default=2)
    ap.add_argument('--top-groups', type=int, default=2)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    queries = read_tsv(args.candidate_queries_tsv)
    positions = load_target_positions(args.target_sites_tsv)
    pos_to_idx = {p: i for i, p in enumerate(positions)}
    scaffold = load_scaffold(args.scaffold_calls_tsv)
    panel = load_panel_candidates(args.panel_candidates_json)

    groups = sorted({rec['group'] for rec in panel if rec['group']})
    if not groups:
        groups = ['UNK']
    # Precompute group-wise counts at every target-site index.
    gpos = {g: [0] * len(positions) for g in groups}
    gneg = {g: [0] * len(positions) for g in groups}
    gall = {g: [0] * len(positions) for g in groups}
    for rec in panel:
        g = rec['group'] if rec['group'] in gpos else groups[0]
        alleles = rec['alleles']
        m = min(len(alleles), len(positions))
        for i in range(m):
            a = alleles[i]
            if a == 1:
                gpos[g][i] += 1; gall[g][i] += 1
            elif a == -1:
                gneg[g][i] += 1; gall[g][i] += 1
    rows = []
    total_support = 0.0
    pos_specific_hits = 0
    top_group_counter = Counter()

    # cache anchor windows by site index to avoid recomputation
    anchor_cache = {}

    for q in queries:
        pos = int(q['pos'])
        cstate = int(q['candidate_phase_state'])
        idx = pos_to_idx.get(pos)
        donor_bias = 0.0
        donor_entropy = 1.0
        candidate_matches_donor = 0
        panel_support_sum = 0.0
        panel_support_frac = 0.0
        panel_oppose_sum = 0.0
        panel_total_informative = 0.0
        panel_top_group = ''
        panel_group_margin = 0.0
        panel_anchor_count = 0

        if idx is not None and cstate != 0:
            if idx not in anchor_cache:
                lo = max(0, idx - args.window_sites)
                hi = min(len(positions) - 1, idx + args.window_sites)
                anchors = []
                for j in range(lo, hi + 1):
                    p = positions[j]
                    if p == pos:
                        continue
                    if p in scaffold:
                        anchors.append((j, scaffold[p][0], scaffold[p][1]))
                anchor_cache[idx] = anchors
            anchors = anchor_cache[idx]
            panel_anchor_count = len(anchors)
            if len(anchors) >= args.min_anchor_count:
                # score each group by local scaffold agreement, weighted by confidence and proximity
                g_scores = {}
                for g in groups:
                    s = 0.0
                    for j, st, conf in anchors:
                        dist = abs(j - idx)
                        w = conf / (1.0 + dist)
                        posc = gpos[g][j]
                        negc = gneg[g][j]
                        if st == 1:
                            s += w * math.log((posc + 1.0) / (negc + 1.0))
                        else:
                            s += w * math.log((negc + 1.0) / (posc + 1.0))
                    g_scores[g] = s
                ranked = sorted(g_scores.items(), key=lambda kv: kv[1], reverse=True)
                selected_groups = [g for g, _ in ranked[:max(1, args.top_groups)]]
                top_group = ranked[0][0]
                top_group_counter[top_group] += 1
                panel_top_group = top_group
                if len(ranked) > 1:
                    panel_group_margin = ranked[0][1] - ranked[1][1]
                probs = softmax([s for _, s in ranked])
                donor_entropy = entropy(probs)

                support = 0.0
                oppose = 0.0
                informative = 0.0
                for g in selected_groups:
                    if cstate == 1:
                        support += gpos[g][idx]
                        oppose += gneg[g][idx]
                    else:
                        support += gneg[g][idx]
                        oppose += gpos[g][idx]
                    informative += gall[g][idx]

                panel_support_sum = float(support)
                panel_oppose_sum = float(oppose)
                panel_total_informative = float(informative)
                if support + oppose > 0:
                    panel_support_frac = support / (support + oppose)
                    donor_bias = (support - oppose) / (support + oppose)
                    candidate_matches_donor = 1 if support > oppose else 0
                    pos_specific_hits += 1
                    total_support += panel_support_sum

        rows.append({
            'query_id': q['query_id'],
            'candidate_id': q['candidate_id'],
            'pos': pos,
            'candidate_phase_state': cstate,
            'donor_bias': f'{donor_bias:.6f}',
            'donor_entropy': f'{donor_entropy:.6f}',
            'candidate_matches_donor': candidate_matches_donor,
            'panel_support_sum': f'{panel_support_sum:.6f}',
            'panel_support_frac': f'{panel_support_frac:.6f}',
            'panel_oppose_sum': f'{panel_oppose_sum:.6f}',
            'panel_total_informative': f'{panel_total_informative:.6f}',
            'panel_top_group': panel_top_group,
            'panel_group_margin': f'{panel_group_margin:.6f}',
            'panel_anchor_count': panel_anchor_count,
        })

    fields = list(rows[0].keys()) if rows else []
    write_tsv(args.out_tsv, rows, fields)
    summary = {
        'rows': len(rows),
        'panel_top_group': top_group_counter.most_common(1)[0][0] if top_group_counter else '',
        'panel_support_sum': total_support,
        'panel_pos_specific_hits': pos_specific_hits,
        'panel_positions': len(positions),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, 'w'), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
