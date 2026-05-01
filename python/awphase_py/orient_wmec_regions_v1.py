#!/usr/bin/env python3
import argparse, csv, json, math
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def to_int(x, default=0):
    try:
        if x is None or str(x).strip() == '':
            return default
        return int(float(x))
    except Exception:
        return default


def to_float(x, default=0.0):
    try:
        if x is None or str(x).strip() == '':
            return default
        return float(x)
    except Exception:
        return default


def load_regions(path):
    regions = {}
    ordered = []
    for r in read_tsv(path):
        rid = r['region_id']
        poss = []
        if r.get('positions'):
            poss = [int(x) for x in str(r['positions']).split(',') if x]
        regions[rid] = {
            'region_id': rid,
            'positions': poss,
            'start': to_int(r.get('start_pos'), min(poss) if poss else 0),
            'end': to_int(r.get('end_pos'), max(poss) if poss else 0),
            'left_anchor_pos': to_int(r.get('left_anchor_pos'), 0),
            'left_anchor_state': to_int(r.get('left_anchor_state'), 0),
            'right_anchor_pos': to_int(r.get('right_anchor_pos'), 0),
            'right_anchor_state': to_int(r.get('right_anchor_state'), 0),
        }
        ordered.append(rid)
    ordered.sort(key=lambda x: (regions[x]['start'], regions[x]['end']))
    return regions, ordered


def load_solutions(path):
    by_region = defaultdict(list)
    pos_to_region = {}
    state_by_pos = {}
    base_rows = []
    for r in read_tsv(path):
        rid = r['region_id']
        pos = to_int(r['pos'])
        st = to_int(r['resolved_state'])
        row = dict(r)
        row['_pos_int'] = pos
        row['_state_int'] = st
        by_region[rid].append(row)
        pos_to_region[pos] = rid
        state_by_pos[pos] = st
        base_rows.append(row)
    for rid in by_region:
        by_region[rid].sort(key=lambda r: r['_pos_int'])
    return by_region, pos_to_region, state_by_pos, base_rows


def load_local_states(path):
    out = {}
    for r in read_tsv(path):
        pos = to_int(r.get('pos'))
        st = to_int(r.get('local_phase_state', r.get('phase_state', 0)))
        if pos:
            out[pos] = st
    return out


def obs_weight(r):
    # Use the best calibrated weight if available; otherwise approximate from confidence/baseQ.
    for k in ('weight_v3','weight_v2','allele_confidence'):
        v = to_float(r.get(k), None)
        if v is not None:
            if k == 'allele_confidence':
                return max(0.05, min(4.0, 1.0 + 6.0 * v))
            return max(0.05, min(6.0, v))
    bq = to_float(r.get('effective_baseq', r.get('raw_baseq', 20)), 20.0)
    return max(0.05, min(6.0, bq / 10.0))


def load_obs_fragments(path, pos_to_region, max_reads=500000):
    # Return read_id -> [(pos, allele, weight, region_id)] for reads touching solved positions.
    reads = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for i, r in enumerate(reader):
            if i > max_reads * 20 and len(reads) > max_reads:
                break
            rid = str(r.get('read_id',''))
            if not rid:
                continue
            pos = to_int(r.get('site_pos', r.get('pos', 0)))
            if pos not in pos_to_region:
                continue
            allele = to_int(r.get('allele', 0))
            if allele not in (-1, 1):
                continue
            if to_int(r.get('is_ambiguous', 0)):
                continue
            w = obs_weight(r)
            reads[rid].append((pos, allele, w, pos_to_region[pos]))
    # keep informative reads only
    return {rid: sorted(v, key=lambda x: x[0]) for rid, v in reads.items() if len({x[3] for x in v}) >= 1}


def fragment_cost(obs, state_by_pos, orient_by_region):
    # Standard fragment-to-haplotype MEC cost: choose the haplotype polarity that minimizes mismatches.
    costs = []
    for hap_pol in (1, -1):
        c = 0.0
        for pos, allele, w, rid in obs:
            oriented_state = state_by_pos.get(pos, 0) * orient_by_region.get(rid, 1)
            if oriented_state == 0:
                continue
            expected = hap_pol * oriented_state
            if allele != expected:
                c += w
        costs.append(c)
    return min(costs) if costs else 0.0


def build_pair_edges(reads, state_by_pos):
    # For every read spanning multiple regions, contribute pairwise orientation costs between adjacent regions on that read.
    edges = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])  # (++,+-,-+,--)
    bridge_reads = defaultdict(int)
    for rid, obs in reads.items():
        regions = []
        by_r = defaultdict(list)
        for o in obs:
            by_r[o[3]].append(o)
        regions = sorted(by_r.keys(), key=lambda r: min(x[0] for x in by_r[r]))
        if len(regions) < 2:
            continue
        for a, b in zip(regions, regions[1:]):
            if a == b:
                continue
            pair_obs = by_r[a] + by_r[b]
            vals = []
            for oa, ob in ((1,1),(1,-1),(-1,1),(-1,-1)):
                vals.append(fragment_cost(pair_obs, state_by_pos, {a: oa, b: ob}))
            key = (a,b)
            for i,v in enumerate(vals):
                edges[key][i] += v
            bridge_reads[key] += 1
    return edges, bridge_reads


def region_unary_costs(rid, region, sol_rows, local_states, anchor_weight, baseline_weight, donor_weight):
    # Cost for orientation +1 and -1. These terms break the inherent MEC orientation symmetry.
    costs = {1: 0.0, -1: 0.0}
    states = [(r['_pos_int'], r['_state_int']) for r in sol_rows]
    if not states:
        return costs

    for orient in (1, -1):
        cost = 0.0
        # Keep close to current selected branch unless evidence says otherwise.
        for pos, st in states:
            local = local_states.get(pos, 0)
            if local and st * orient != local:
                cost += baseline_weight

        # Anchor terms: lightweight but important for phase polarity.
        if region.get('left_anchor_state'):
            if states[0][1] * orient != region['left_anchor_state']:
                cost += anchor_weight
        if region.get('right_anchor_state'):
            if states[-1][1] * orient != region['right_anchor_state']:
                cost += anchor_weight

        # If original solve says anchor_ok, flipping should require evidence.
        anchor_ok = sum(to_int(r.get('anchor_ok', 0)) for r in sol_rows)
        if orient == -1 and anchor_ok > 0:
            cost += anchor_weight * min(2.0, anchor_ok / max(1, len(sol_rows)))

        # Donor term: avg_donor_bias is directional for the reported resolved state.
        avg_donor = sum(to_float(r.get('avg_donor_bias', 0.0)) for r in sol_rows) / max(1, len(sol_rows))
        cost += - donor_weight * orient * avg_donor
        costs[orient] = cost
    return costs


def chain_dp(active_regions, unaries, edges):
    # DP over genomic order; transitions only for adjacent active regions, with pair edges if present.
    if not active_regions:
        return {}, {}
    dp = []
    back = []
    first = active_regions[0]
    dp.append({1: unaries[first][1], -1: unaries[first][-1]})
    back.append({1: None, -1: None})
    for i in range(1, len(active_regions)):
        prev = active_regions[i-1]
        cur = active_regions[i]
        ec = edges.get((prev, cur), None)
        curdp = {}
        curback = {}
        for oc in (1, -1):
            best = None
            bestop = None
            for op in (1, -1):
                trans = 0.0
                if ec is not None:
                    idx = {(1,1):0,(1,-1):1,(-1,1):2,(-1,-1):3}[(op,oc)]
                    trans = ec[idx]
                val = dp[-1][op] + trans + unaries[cur][oc]
                if best is None or val < best:
                    best = val; bestop = op
            curdp[oc] = best; curback[oc] = bestop
        dp.append(curdp); back.append(curback)

    last_orient = 1 if dp[-1][1] <= dp[-1][-1] else -1
    best_cost = dp[-1][last_orient]
    orient = {}
    cur_o = last_orient
    for i in range(len(active_regions)-1, -1, -1):
        rid = active_regions[i]
        orient[rid] = cur_o
        cur_o = back[i][cur_o]
        if cur_o is None:
            break
    orient = {rid: orient.get(rid, 1) for rid in active_regions}

    # Compute per-region local orientation margin from unary + adjacent edge approximations.
    margins = {}
    for rid in active_regions:
        c1 = unaries[rid][1]
        c2 = unaries[rid][-1]
        margins[rid] = abs(c1 - c2)
    return orient, margins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--regions-tsv', required=True)
    ap.add_argument('--wmec-solutions-tsv', required=True)
    ap.add_argument('--local-calls-tsv', required=True)
    ap.add_argument('--obs-tsv', required=True)
    ap.add_argument('--out-oriented-tsv', required=True)
    ap.add_argument('--out-orientations-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    ap.add_argument('--anchor-weight', type=float, default=2.5)
    ap.add_argument('--baseline-weight', type=float, default=0.25)
    ap.add_argument('--donor-weight', type=float, default=1.0)
    ap.add_argument('--min-flip-margin', type=float, default=0.75)
    args = ap.parse_args()

    regions, ordered_all = load_regions(args.regions_tsv)
    sol_by_region, pos_to_region, state_by_pos, sol_rows = load_solutions(args.wmec_solutions_tsv)
    local_states = load_local_states(args.local_calls_tsv)
    reads = load_obs_fragments(args.obs_tsv, pos_to_region)
    edges, bridge_reads = build_pair_edges(reads, state_by_pos)

    active_regions = [rid for rid in ordered_all if rid in sol_by_region]
    unaries = {}
    for rid in active_regions:
        unaries[rid] = region_unary_costs(rid, regions.get(rid, {}), sol_by_region[rid], local_states, args.anchor_weight, args.baseline_weight, args.donor_weight)

    orient_raw, margins = chain_dp(active_regions, unaries, edges)

    # Conservative orientation: only allow a flip if the local unary margin is meaningful or there is bridge support.
    final_orient = {}
    for rid in active_regions:
        o = orient_raw.get(rid, 1)
        has_bridge = any(rid in pair for pair in edges.keys())
        if o == -1 and margins.get(rid, 0.0) < args.min_flip_margin and not has_bridge:
            o = 1
        final_orient[rid] = o

    orient_rows = []
    for rid in active_regions:
        pair_support = sum(v for pair,v in bridge_reads.items() if rid in pair)
        orient_rows.append({
            'region_id': rid,
            'orientation': final_orient.get(rid,1),
            'raw_orientation': orient_raw.get(rid,1),
            'orientation_margin': f"{margins.get(rid,0.0):.6f}",
            'unary_plus': f"{unaries[rid][1]:.6f}",
            'unary_minus': f"{unaries[rid][-1]:.6f}",
            'bridge_read_pairs': pair_support,
            'n_sites': len(sol_by_region[rid]),
        })

    oriented_rows = []
    for r in sol_rows:
        row = {k:v for k,v in r.items() if not k.startswith('_')}
        rid = row['region_id']
        o = final_orient.get(rid, 1)
        old = to_int(row.get('resolved_state', 0))
        row['resolved_state_original'] = old
        row['resolved_state'] = old * o
        row['region_orientation'] = o
        row['orientation_margin'] = f"{margins.get(rid,0.0):.6f}"
        oriented_rows.append(row)

    write_tsv(args.out_orientations_tsv, orient_rows, list(orient_rows[0].keys()) if orient_rows else ['region_id'])
    write_tsv(args.out_oriented_tsv, oriented_rows, list(oriented_rows[0].keys()) if oriented_rows else ['region_id','pos','resolved_state'])

    summary = {
        'regions': len(active_regions),
        'flipped_regions': sum(1 for r in active_regions if final_orient.get(r,1) == -1),
        'raw_flipped_regions': sum(1 for r in active_regions if orient_raw.get(r,1) == -1),
        'bridge_edges': len(edges),
        'reads_touching_solved_positions': len(reads),
        'oriented_rows': len(oriented_rows),
        'min_flip_margin': args.min_flip_margin,
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json,'w'), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
