#!/usr/bin/env python3
import argparse, csv, json
from collections import defaultdict, deque, Counter
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def as_int(x, d=0):
    try: return int(float(x))
    except Exception: return d


def as_float(x, d=0.0):
    try: return float(x)
    except Exception: return d


def detect_score_col(row):
    for c in ("rank_score","pred_score","score","candidate_score","xgb_score","prob","ranker_score"):
        if c in row: return c
    for c in row:
        if "score" in c.lower(): return c
    raise SystemExit(f"Could not find score column in candidate scores: {list(row.keys())}")


def load_candidate_scores(path):
    rows = read_tsv(path)
    if not rows:
        raise SystemExit(f"Empty candidate scores TSV: {path}")
    score_col = detect_score_col(rows[0])
    by_pos = defaultdict(list)
    for r in rows:
        pos = as_int(r.get("pos", r.get("site_pos")))
        if not pos: continue
        state = as_int(r.get("candidate_phase_state", r.get("phase_state")))
        by_pos[pos].append({
            "candidate_id": r.get("candidate_id",""),
            "state": state,
            "score": as_float(r.get(score_col,0.0)),
            "candidate_type": r.get("candidate_type",""),
            "donor_bias": as_float(r.get("f_donor_bias", r.get("donor_bias", 0.0))),
            "anchor_agreement": as_float(r.get("f_anchor_agreement", r.get("anchor_agreement", 0.0))),
        })
    pos_info = {}
    for pos, lst in by_pos.items():
        lst = sorted(lst, key=lambda x: x["score"], reverse=True)
        best = lst[0]
        second = lst[1]["score"] if len(lst) > 1 else -999.0
        non_abstain = [x for x in lst if x["state"] != 0]
        best_non = non_abstain[0] if non_abstain else best
        pos_info[pos] = {
            "top_state": best["state"],
            "top_nonzero_state": best_non["state"],
            "top_score": best["score"],
            "margin": best["score"] - second,
            "top_donor_bias": best["donor_bias"],
            "top_anchor_agreement": best["anchor_agreement"],
            "has_non_abstain": bool(non_abstain),
        }
    return pos_info


def load_local(path):
    out = {}
    rows = read_tsv(path)
    if not rows: return out
    pos_col = "pos" if "pos" in rows[0] else next(iter(rows[0].keys()))
    state_col = "local_phase_state" if "local_phase_state" in rows[0] else "phase_state"
    for r in rows:
        pos = as_int(r.get(pos_col))
        if not pos: continue
        out[pos] = {
            "local_state": as_int(r.get(state_col, 0)),
            "confidence": as_float(r.get("confidence", r.get("local_confidence", 0.0))),
            "block_id": r.get("block_id", r.get("scaffold_block_id","")),
        }
    return out


def load_scaffold(path):
    rows = read_tsv(path)
    out = {}; scaffold_positions = []
    for r in rows:
        pos = as_int(r.get("pos", r.get("site_pos")))
        if not pos: continue
        state = as_int(r.get("scaffold_phase_state", 0))
        conf = as_float(r.get("scaffold_confidence", 0.0))
        bid = r.get("scaffold_block_id","")
        out[pos] = {"state": state, "conf": conf, "block_id": bid}
        if state != 0:
            scaffold_positions.append(pos)
    scaffold_positions.sort()
    return out, scaffold_positions


def load_obs_graph(path, min_conf, min_delta, min_sites_per_read, max_link_gap_bp):
    # Build full read-connectivity graph, not only candidate-to-candidate adjacency.
    read_to_positions = defaultdict(set)
    connected_reads = 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            pos = as_int(r.get("site_pos", r.get("pos")))
            allele = as_int(r.get("allele", 0))
            ambig = as_int(r.get("is_ambiguous", 0))
            if not pos or allele not in (-1, 1) or ambig != 0:
                continue
            conf = as_float(r.get("allele_confidence", 0.0))
            delta = abs(as_float(r.get("allele_score_delta", 0.0)))
            # Cleaner matrix but not starved: accept if either signal is strong enough.
            if conf < min_conf and delta < min_delta:
                continue
            rid = r.get("read_id", "").strip()
            if rid:
                read_to_positions[rid].add(pos)
    adj = defaultdict(set)
    read_links = {}
    for rid, posset in read_to_positions.items():
        ps = sorted(posset)
        if len(ps) < min_sites_per_read:
            continue
        connected_reads += 1
        read_links[rid] = ps
        # Link all consecutive variants spanned by a read, but avoid absurd long joins.
        for a, b in zip(ps, ps[1:]):
            if (b - a) <= max_link_gap_bp:
                adj[a].add(b); adj[b].add(a)
    return adj, read_links, connected_reads


def nearest_anchor(scaffold_positions, scaffold_map, pos, side):
    if side == "left":
        lo = 0; hi = len(scaffold_positions)
        # simple binary search by hand to avoid importing bisect? no, use loop OK at this scale.
        cands = [p for p in scaffold_positions if p < pos]
        if not cands: return ("", 0)
        p = cands[-1]
    else:
        cands = [p for p in scaffold_positions if p > pos]
        if not cands: return ("", 0)
        p = cands[0]
    return (p, scaffold_map[p]["state"])


def split_component(component, max_sites, max_span):
    comp = sorted(component)
    chunks = []
    cur = []; start = None
    for p in comp:
        if not cur:
            cur = [p]; start = p; continue
        if len(cur) >= max_sites or (p - start) > max_span:
            chunks.append(cur); cur = [p]; start = p
        else:
            cur.append(p)
    if cur: chunks.append(cur)
    return chunks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--candidate-scores-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--scaffold-calls-tsv", required=True)
    ap.add_argument("--min-candidate-margin", type=float, default=0.02)
    ap.add_argument("--min-obs-confidence", type=float, default=0.025)
    ap.add_argument("--min-obs-delta", type=float, default=0.45)
    ap.add_argument("--min-sites-per-read", type=int, default=2)
    ap.add_argument("--max-sites-per-region", type=int, default=24)
    ap.add_argument("--max-region-span-bp", type=int, default=100000)
    ap.add_argument("--max-read-link-gap-bp", type=int, default=100000)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    pos_info = load_candidate_scores(args.candidate_scores_tsv)
    local = load_local(args.local_calls_tsv)
    scaffold, scaffold_positions = load_scaffold(args.scaffold_calls_tsv)
    adj, read_links, connected_reads = load_obs_graph(args.obs_tsv, args.min_obs_confidence, args.min_obs_delta, args.min_sites_per_read, args.max_read_link_gap_bp)

    candidate_positions = set()
    for pos, info in pos_info.items():
        local_state = local.get(pos, {}).get("local_state", 0)
        top = info["top_state"]
        # Candidate if ranker proposes a non-local state, if local abstained, or if confidence margin is meaningful.
        if (top != 0 and top != local_state) or local_state == 0 or info["margin"] >= args.min_candidate_margin:
            candidate_positions.add(pos)

    # Walk full graph components and collect candidate positions reached through non-candidate bridge nodes.
    graph_nodes = set(adj.keys()) | {v for vs in adj.values() for v in vs}
    visited = set(); regions = []
    for node in sorted(graph_nodes):
        if node in visited: continue
        dq = deque([node]); visited.add(node); comp_nodes = []
        while dq:
            u = dq.popleft(); comp_nodes.append(u)
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v); dq.append(v)
        cand = sorted(set(comp_nodes).intersection(candidate_positions))
        if cand:
            regions.extend(split_component(cand, args.max_sites_per_region, args.max_region_span_bp))

    # Add candidate positions absent from the graph as singleton regions.
    assigned = {p for reg in regions for p in reg}
    for p in sorted(candidate_positions - assigned):
        regions.append([p])
    regions.sort(key=lambda x: (x[0], x[-1], len(x)))

    rows = []
    size_counter = Counter()
    for idx, region in enumerate(regions, start=1):
        left_pos, left_state = nearest_anchor(scaffold_positions, scaffold, region[0], "left")
        right_pos, right_state = nearest_anchor(scaffold_positions, scaffold, region[-1], "right")
        size_counter[len(region)] += 1
        rows.append({
            "region_id": f"rg_{idx}",
            "n_sites": len(region),
            "start_pos": region[0],
            "end_pos": region[-1],
            "span_bp": region[-1] - region[0] if len(region) > 1 else 0,
            "positions": ",".join(str(p) for p in region),
            "left_anchor_pos": left_pos,
            "left_anchor_state": left_state,
            "right_anchor_pos": right_pos,
            "right_anchor_state": right_state,
        })

    fields = ["region_id","n_sites","start_pos","end_pos","span_bp","positions","left_anchor_pos","left_anchor_state","right_anchor_pos","right_anchor_state"]
    write_tsv(args.out_tsv, rows, fields)
    summary = {
        "regions": len(rows),
        "candidate_positions": len(candidate_positions),
        "connected_reads": connected_reads,
        "graph_nodes": len(graph_nodes),
        "min_candidate_margin": args.min_candidate_margin,
        "max_sites_per_region": args.max_sites_per_region,
        "singleton_regions": size_counter.get(1,0),
        "multi_site_regions": sum(v for k,v in size_counter.items() if k > 1),
        "region_size_hist": dict(sorted(size_counter.items())),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
