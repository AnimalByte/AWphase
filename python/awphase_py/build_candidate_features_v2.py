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

def f(x, d=0.0):
    try: return float(x)
    except Exception: return d

def i(x, d=0):
    try: return int(float(x))
    except Exception: return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-features-tsv', required=True)
    ap.add_argument('--donor-prior-v2-tsv', required=True)
    ap.add_argument('--realign-features-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    base = read_tsv(args.base_features_tsv)
    donor = {(r['query_id'], r['candidate_id']): r for r in read_tsv(args.donor_prior_v2_tsv)}
    realign = {(r['query_id'], r['candidate_id']): r for r in read_tsv(args.realign_features_tsv)}

    scaffold = {}
    for r in read_tsv(args.scaffold_calls_tsv):
        try:
            pos = int(r['pos'])
            sp = int(float(r.get('scaffold_phase_state', 0) or 0))
            conf = float(r.get('scaffold_confidence', 0.0) or 0.0)
        except Exception:
            continue
        scaffold[pos] = (sp, conf)

    scaffold_positions = sorted(scaffold)

    def nearest_left_right(pos):
        left = None; right = None
        for p in scaffold_positions:
            if p < pos: left = p
            elif p > pos:
                right = p
                break
        return left, right

    out = []
    for r in base:
        pos = int(r['pos']); cstate = int(r['candidate_phase_state'])
        key = (r['query_id'], r['candidate_id'])
        d = donor.get(key, {})
        rl = realign.get(key, {})

        left, right = nearest_left_right(pos)
        lstate = scaffold.get(left, (0,0.0))[0] if left is not None else 0
        rstate = scaffold.get(right, (0,0.0))[0] if right is not None else 0
        lconf = scaffold.get(left, (0,0.0))[1] if left is not None else 0.0
        rconf = scaffold.get(right, (0,0.0))[1] if right is not None else 0.0

        matches_left = 1 if (left is not None and cstate != 0 and cstate == lstate) else 0
        matches_right = 1 if (right is not None and cstate != 0 and cstate == rstate) else 0
        anchor_agree = 1 if (left is not None and right is not None and cstate != 0 and cstate == lstate == rstate) else 0
        anchor_conflict = 1 if (left is not None and right is not None and lstate != 0 and rstate != 0 and lstate != rstate) else 0

        row = dict(r)
        row.update({
            'f_donor_bias_local': d.get('donor_bias_local', '0.0'),
            'f_donor_entropy_local': d.get('donor_entropy_local', '1.0'),
            'f_candidate_matches_donor_local': d.get('candidate_matches_donor_local', '0'),
            'f_panel_support_sum_local': d.get('panel_support_sum_local', '0.0'),
            'f_panel_support_frac_local': d.get('panel_support_frac_local', '0.0'),
            'f_panel_oppose_sum_local': d.get('panel_oppose_sum_local', '0.0'),
            'f_panel_total_informative_local': d.get('panel_total_informative_local', '0.0'),
            'f_top_donor_mean_score_local': d.get('top_donor_mean_score_local', '0.0'),
            'f_donor_window_anchor_support': d.get('donor_window_anchor_support', '0.0'),
            'f_donor_window_anchor_disagree': d.get('donor_window_anchor_disagree', '0.0'),
            'f_donor_window_entropy': d.get('donor_window_entropy', '1.0'),
            'f_realign_support_w': rl.get('realign_support_w', '0.0'),
            'f_realign_oppose_w': rl.get('realign_oppose_w', '0.0'),
            'f_realign_margin': rl.get('realign_margin', '0.0'),
            'f_realign_support_delta_sum': rl.get('realign_support_delta_sum', '0.0'),
            'f_realign_oppose_delta_sum': rl.get('realign_oppose_delta_sum', '0.0'),
            'f_realign_tagged_obs': rl.get('realign_tagged_obs', '0'),
            'f_matches_left_anchor': matches_left,
            'f_matches_right_anchor': matches_right,
            'f_anchor_consistent': anchor_agree,
            'f_anchor_conflict': anchor_conflict,
            'f_left_anchor_conf': f'{lconf:.6f}',
            'f_right_anchor_conf': f'{rconf:.6f}',
            'f_local_vs_anchor_disagree': 1 if cstate != 0 and ((left is not None and lstate != 0 and cstate != lstate) or (right is not None and rstate != 0 and cstate != rstate)) else 0,
            'f_nonabstain': 0 if cstate == 0 else 1,
        })
        out.append(row)

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else [])
    summary = {'rows': len(out), 'queries': len({r['query_id'] for r in out})}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
