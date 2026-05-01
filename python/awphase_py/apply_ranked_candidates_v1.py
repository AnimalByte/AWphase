#!/usr/bin/env python3
import argparse, csv
from pathlib import Path
from collections import defaultdict


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    base_fields = list(fields) if fields is not None else []
    seen = set(base_fields)

    # append any new keys introduced during application
    for r in rows:
        for k in r.keys():
            if k not in seen:
                base_fields.append(k)
                seen.add(k)

    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=base_fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def as_int(x, d=0):
    try: return int(float(x))
    except Exception: return d


def as_float(x, d=0.0):
    try: return float(x)
    except Exception: return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--local-calls-tsv', required=True)
    ap.add_argument('--scored-candidates-tsv', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--min-best-score', type=float, default=0.60)
    ap.add_argument('--min-margin', type=float, default=0.05)
    ap.add_argument('--min-modify-score', type=float, default=0.72)
    ap.add_argument('--min-modify-margin', type=float, default=0.15)
    ap.add_argument('--anchor-agreement-min', type=float, default=0.0)
    ap.add_argument('--strong-isolated-score', type=float, default=0.85)
    ap.add_argument('--strong-isolated-margin', type=float, default=0.25)
    ap.add_argument('--max-gap-bp', type=int, default=5000)
    args = ap.parse_args()

    local = read_tsv(args.local_calls_tsv)
    scored = read_tsv(args.scored_candidates_tsv)
    by_q = defaultdict(list)
    for r in scored:
        by_q[r['query_id']].append(r)

    proposals = []
    proposal_by_pos = {}
    for qid, rows in by_q.items():
        rows = sorted(rows, key=lambda r: as_float(r.get('pred_score', 0.0)), reverse=True)
        top = rows[0]
        best_score = as_float(top.get('pred_score', 0.0))
        margin = as_float(top.get('query_margin', 0.0))
        cand_state = as_int(top.get('candidate_phase_state', 0))
        pos = as_int(top['pos'])
        local_state = as_int(top.get('f_candidate_matches_local', 0))
        current_state = None
        # recover current local state from candidates using matches_local signal
        for r in rows:
            if as_int(r.get('f_candidate_matches_local',0)) == 1:
                current_state = as_int(r.get('candidate_phase_state',0))
                break
        if current_state is None:
            current_state = 0
        modify = cand_state != current_state
        anchor_ag = as_float(top.get('f_anchor_agreement', 0.0))
        switch_risk = as_int(top.get('f_anchor_window_switch_risk', 0))
        donor_bias = as_float(top.get('f_donor_bias', 0.0))
        donor_entropy = as_float(top.get('f_donor_entropy', 1.0))
        nearest_anchor = as_float(top.get('f_nearest_anchor_dist_bp', 1e9))

        accept = True
        reason = 'accept'
        if cand_state == 0:
            accept = False; reason = 'top_is_abstain'
        elif best_score < args.min_best_score or margin < args.min_margin:
            accept = False; reason = 'weak_top'
        elif modify and (best_score < args.min_modify_score or margin < args.min_modify_margin):
            accept = False; reason = 'weak_modify'
        elif modify and (anchor_ag < args.anchor_agreement_min) and nearest_anchor < 1e9:
            accept = False; reason = 'anchor_disagree'
        elif modify and switch_risk == 1 and margin < max(args.min_modify_margin, 0.20):
            accept = False; reason = 'switch_risk_guard'
        elif modify and donor_bias < -0.05:
            accept = False; reason = 'donor_disagree'
        elif modify and donor_entropy > 0.95 and margin < 0.20:
            accept = False; reason = 'donor_uncertain'

        if accept:
            proposal = {
                'query_id': qid,
                'pos': pos,
                'current_state': current_state,
                'proposed_state': cand_state,
                'best_score': best_score,
                'margin': margin,
                'anchor_agreement': anchor_ag,
                'donor_bias': donor_bias,
                'reason': reason,
            }
            proposals.append(proposal)
            proposal_by_pos[pos] = proposal

    # suppress isolated flips unless very strong
    proposals.sort(key=lambda x: x['pos'])
    keep_positions = set()
    i = 0
    while i < len(proposals):
        j = i + 1
        run = [proposals[i]]
        while j < len(proposals) and proposals[j]['proposed_state'] == run[-1]['proposed_state'] and proposals[j]['pos'] - run[-1]['pos'] <= args.max_gap_bp:
            run.append(proposals[j]); j += 1
        if len(run) == 1:
            p = run[0]
            if p['best_score'] >= args.strong_isolated_score and p['margin'] >= args.strong_isolated_margin and (p['anchor_agreement'] > 0 or p['donor_bias'] > 0.10):
                keep_positions.add(p['pos'])
        else:
            avg_margin = sum(r['margin'] for r in run) / len(run)
            if avg_margin >= args.min_modify_margin:
                for r in run:
                    keep_positions.add(r['pos'])
        i = j

    out = []
    modified = 0
    for r in local:
        row = dict(r)
        pos = as_int(r['pos'])
        prop = proposal_by_pos.get(pos)
        if prop and pos in keep_positions:
            old = as_int(r.get('local_phase_state', 0))
            new = prop['proposed_state']
            if new != old:
                row['local_phase_state'] = str(new)
                modified += 1
                row['ranked_applied'] = '1'
                row['ranked_best_score'] = f"{prop['best_score']:.6f}"
                row['ranked_margin'] = f"{prop['margin']:.6f}"
                row['ranked_anchor_agreement'] = f"{prop['anchor_agreement']:.6f}"
                row['ranked_donor_bias'] = f"{prop['donor_bias']:.6f}"
                row['ranked_reason'] = prop['reason']
            else:
                row['ranked_applied'] = '0'
        else:
            row['ranked_applied'] = '0'
        out.append(row)

    fields = list(out[0].keys()) if out else []
    write_tsv(args.out_tsv, out, fields)
    print({'rows': len(out), 'modified_rows': modified, 'threshold': args.min_best_score})

if __name__ == '__main__':
    main()
