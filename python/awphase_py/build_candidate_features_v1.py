#!/usr/bin/env python3
import argparse, csv, json
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)


def as_int(x, d=0):
    try: return int(float(x))
    except Exception: return d


def as_float(x, d=0.0):
    try: return float(x)
    except Exception: return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-queries-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--donor-prior-tsv', required=True)
    ap.add_argument('--reference-context-tsv', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    qrows = read_tsv(args.candidate_queries_tsv)
    srows = read_tsv(args.scaffold_calls_tsv)
    drows = read_tsv(args.donor_prior_tsv)
    rrows = read_tsv(args.reference_context_tsv)

    scaffold = {}
    anchors = []
    for r in srows:
        pos = as_int(r['pos'])
        st = as_int(r.get('scaffold_phase_state', 0))
        record = {
            'local_phase_state': as_int(r.get('local_phase_state', 0)),
            'local_confidence': as_float(r.get('local_confidence', 0.0)),
            'scaffold_phase_state': st,
            'scaffold_confidence': as_float(r.get('scaffold_confidence', 0.0)),
            'selected_support': as_float(r.get('selected_support', 0.0)),
            'tagged_support': as_float(r.get('tagged_support', 0.0)),
            'superread_support': as_float(r.get('superread_support', 0.0)),
            'site_prior_bias': as_float(r.get('site_prior_bias', 0.0)),
            'ambiguous_frac': as_float(r.get('ambiguous_frac', 0.0)),
            'mean_allele_confidence': as_float(r.get('mean_allele_confidence', 0.0)),
            'mean_abs_score_delta': as_float(r.get('mean_abs_score_delta', 0.0)),
            'hp1_reads': as_float(r.get('hp1_reads', 0.0)),
            'hp2_reads': as_float(r.get('hp2_reads', 0.0)),
            'scaffold_block_id': r.get('scaffold_block_id', ''),
        }
        scaffold[pos] = record
        if st != 0:
            anchors.append((pos, st, record['scaffold_confidence']))
    anchors.sort()

    donor = {(r['query_id'], r['candidate_id']): r for r in drows}
    refctx = {as_int(r['pos']): r for r in rrows}

    # precompute nearest anchors for all positions seen in queries
    q_positions = sorted({as_int(r['pos']) for r in qrows})
    left_map, right_map = {}, {}
    ai = 0
    last = None
    for p in q_positions:
        while ai < len(anchors) and anchors[ai][0] < p:
            last = anchors[ai]
            ai += 1
        left_map[p] = last
    ai = len(anchors) - 1
    nxt = None
    for p in reversed(q_positions):
        while ai >= 0 and anchors[ai][0] > p:
            nxt = anchors[ai]
            ai -= 1
        right_map[p] = nxt

    out = []
    for q in qrows:
        pos = as_int(q['pos'])
        cstate = as_int(q.get('candidate_phase_state', 0))
        s = scaffold.get(pos, {})
        d = donor.get((q['query_id'], q['candidate_id']), {})
        rc = refctx.get(pos, {})

        left = left_map.get(pos)
        right = right_map.get(pos)
        has_left = 1 if left else 0
        has_right = 1 if right else 0
        left_dist = (pos - left[0]) if left else 1e9
        right_dist = (right[0] - pos) if right else 1e9
        left_match = 1 if left and cstate != 0 and left[1] == cstate else 0
        right_match = 1 if right and cstate != 0 and right[1] == cstate else 0
        left_dis = 1 if left and cstate != 0 and left[1] != cstate else 0
        right_dis = 1 if right and cstate != 0 and right[1] != cstate else 0
        n_anchor = has_left + has_right
        anchor_agreement = ((left_match + right_match) - (left_dis + right_dis)) / n_anchor if n_anchor else 0.0
        anchor_switch_risk = 1 if left and right and left[1] != right[1] else 0
        nearest_anchor_dist = min(left_dist, right_dist)
        total_read_support = s.get('selected_support', 0.0) + s.get('tagged_support', 0.0) + s.get('superread_support', 0.0)
        hp_total = s.get('hp1_reads',0.0) + s.get('hp2_reads',0.0)
        hp_imbalance = abs(s.get('hp1_reads',0.0) - s.get('hp2_reads',0.0)) / hp_total if hp_total else 0.0
        candidate_matches_local = 1 if cstate != 0 and cstate == s.get('local_phase_state', 0) else 0

        row = {
            'query_id': q['query_id'],
            'candidate_id': q['candidate_id'],
            'pos': pos,
            'candidate_phase_state': cstate,
            'candidate_type': q.get('candidate_type', ''),
            'f_local_confidence': f"{s.get('local_confidence', 0.0):.6f}",
            'f_scaffold_confidence': f"{s.get('scaffold_confidence', 0.0):.6f}",
            'f_selected_support': f"{s.get('selected_support', 0.0):.6f}",
            'f_tagged_support': f"{s.get('tagged_support', 0.0):.6f}",
            'f_superread_support': f"{s.get('superread_support', 0.0):.6f}",
            'f_total_read_support': f"{total_read_support:.6f}",
            'f_hp_imbalance': f"{hp_imbalance:.6f}",
            'f_site_prior_bias': f"{s.get('site_prior_bias', 0.0):.6f}",
            'f_abs_site_prior_bias': f"{abs(s.get('site_prior_bias', 0.0)):.6f}",
            'f_ambiguous_frac': f"{s.get('ambiguous_frac', 0.0):.6f}",
            'f_mean_allele_confidence': f"{s.get('mean_allele_confidence', 0.0):.6f}",
            'f_mean_abs_score_delta': f"{s.get('mean_abs_score_delta', 0.0):.6f}",
            'f_donor_bias': d.get('donor_bias', '0.0'),
            'f_donor_entropy': d.get('donor_entropy', '1.0'),
            'f_candidate_matches_donor': d.get('candidate_matches_donor', '0'),
            'f_abstain_bonus': '0.500000' if cstate == 0 else '0.000000',
            'f_panel_support_sum': d.get('panel_support_sum', '0.0'),
            'f_panel_support_frac': d.get('panel_support_frac', '0.0'),
            'f_panel_oppose_sum': d.get('panel_oppose_sum', '0.0'),
            'f_panel_total_informative': d.get('panel_total_informative', '0.0'),
            'f_panel_group_margin': d.get('panel_group_margin', '0.0'),
            'f_panel_anchor_count': d.get('panel_anchor_count', '0'),
            'f_anchor_agreement': f"{anchor_agreement:.6f}",
            'f_anchor_window_switch_risk': anchor_switch_risk,
            'f_left_anchor_dist_bp': left_dist,
            'f_right_anchor_dist_bp': right_dist,
            'f_nearest_anchor_dist_bp': nearest_anchor_dist,
            'f_has_left_anchor': has_left,
            'f_has_right_anchor': has_right,
            'f_matches_left_anchor': left_match,
            'f_matches_right_anchor': right_match,
            'f_candidate_matches_local': candidate_matches_local,
            'f_candidate_is_abstain': 1 if cstate == 0 else 0,
            'f_candidate_is_h1': 1 if cstate == 1 else 0,
            'f_candidate_is_h2': 1 if cstate == -1 else 0,
            'f_ctx_gc_frac': rc.get('ctx_gc_frac', '0.0'),
            'f_ctx_entropy': rc.get('ctx_entropy', '0.0'),
            'f_ctx_longest_hp': rc.get('ctx_longest_hp', '0'),
            'f_left_hp': rc.get('left_hp', '0'),
            'f_right_hp': rc.get('right_hp', '0'),
            'f_neighbor_dist_prev': rc.get('neighbor_dist_prev', '1000000000'),
            'f_neighbor_dist_next': rc.get('neighbor_dist_next', '1000000000'),
            'f_low_complexity_flag': rc.get('low_complexity_flag', '0'),
            'f_is_indel': rc.get('is_indel', '0'),
        }
        out.append(row)

    fields = list(out[0].keys()) if out else ['query_id']
    write_tsv(args.out_tsv, out, fields)
    summ = {'rows': len(out), 'queries': len({r['query_id'] for r in out}), 'reference_rows': len(refctx)}
    json.dump(summ, open(args.out_summary_json, 'w'), indent=2)
    print(json.dumps(summ, indent=2))

if __name__ == '__main__':
    main()
