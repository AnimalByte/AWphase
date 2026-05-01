#!/usr/bin/env python3
import argparse, csv, json, math
from collections import defaultdict, Counter
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


def parse_positions(s):
    if not s: return []
    return [as_int(x) for x in str(s).split(',') if str(x).strip()]


def load_regions(path):
    out = []
    for r in read_tsv(path):
        ps = parse_positions(r.get("positions", ""))
        if not ps: continue
        out.append({
            "region_id": r["region_id"],
            "positions": ps,
            "left_anchor_state": as_int(r.get("left_anchor_state", 0)),
            "right_anchor_state": as_int(r.get("right_anchor_state", 0)),
        })
    return out


def load_obs(obs_tsv, min_conf, min_delta):
    read_obs = defaultdict(dict)
    with open(obs_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            rid = r.get("read_id", "").strip()
            pos = as_int(r.get("site_pos", r.get("pos")))
            allele = as_int(r.get("allele", 0))
            if not rid or not pos or allele not in (-1, 1):
                continue
            ambig = as_int(r.get("is_ambiguous", 0))
            if ambig != 0:
                continue
            conf = as_float(r.get("allele_confidence", 0.0))
            delta = abs(as_float(r.get("allele_score_delta", 0.0)))
            if conf < min_conf and delta < min_delta:
                continue
            basew = as_float(r.get("weight_v3", r.get("weight_v2", 0.5)))
            effbq = as_float(r.get("effective_baseq", 20.0)) / 40.0
            w = max(0.15, min(8.0, 1.5*basew + 2.0*conf + min(delta, 3.0) + effbq))
            prev = read_obs[rid].get(pos)
            if prev is None or prev[1] < w:
                read_obs[rid][pos] = (allele, w)
    return read_obs


def fragment_score(posmap):
    ps = sorted(posmap)
    if not ps: return 0.0
    span = ps[-1] - ps[0] + 1
    avgw = sum(w for _, w in posmap.values()) / len(posmap)
    return len(ps) * avgw * math.log(span + 10.0)


def cap_reads_for_region(region_positions, read_obs, max_cov, max_reads):
    rset = set(region_positions)
    candidate = []
    # Allow single-site fragments only for singleton regions; otherwise require bridge reads.
    min_sites = 1 if len(region_positions) == 1 else 2
    for rid, obs in read_obs.items():
        sub = {p: v for p, v in obs.items() if p in rset}
        if len(sub) >= min_sites:
            candidate.append((fragment_score(sub), rid, sub))
    candidate.sort(reverse=True)
    cov = defaultdict(int)
    selected = []
    for _, rid, sub in candidate:
        if len(selected) >= max_reads:
            break
        # Add read if it contributes to at least one under-covered site.
        if any(cov[p] < max_cov for p in sub):
            selected.append((rid, sub))
            for p in sub:
                cov[p] += 1
    return selected


def same_haplotype_logodds(a, b, p_same=0.97, p_diff=0.03):
    overlap = set(a).intersection(b)
    if len(overlap) < 2:
        return -999.0
    agree = sum(1 for p in overlap if a[p][0] == b[p][0])
    disagree = len(overlap) - agree
    ll_same = agree * math.log(p_same) + disagree * math.log(1.0 - p_same)
    ll_diff = agree * math.log(p_diff) + disagree * math.log(1.0 - p_diff)
    return ll_same - ll_diff


def merge_probabilistically(selected, merge_logodds):
    # Stricter than v1: merge only strong same-haplotype evidence. This preserves coverage diversity.
    clusters = []
    for rid, sub in selected:
        best_idx = -1; best_lo = merge_logodds
        for i, cl in enumerate(clusters):
            lo = same_haplotype_logodds(sub, cl["posmap"])
            if lo > best_lo:
                best_lo = lo; best_idx = i
        if best_idx == -1:
            clusters.append({"members": [rid], "posmap": dict(sub)})
        else:
            cl = clusters[best_idx]
            cl["members"].append(rid)
            for p, (a, w) in sub.items():
                if p not in cl["posmap"]:
                    cl["posmap"][p] = (a, w)
                else:
                    a0, w0 = cl["posmap"][p]
                    if a0 == a:
                        cl["posmap"][p] = (a0, w0 + w)
                    else:
                        # weighted cancellation; retain uncertainty but don't flip on small evidence
                        if w > w0:
                            cl["posmap"][p] = (a, max(0.01, w - w0))
                        else:
                            cl["posmap"][p] = (a0, max(0.01, w0 - w))
    return clusters


def load_scaffold_calls(path):
    rows = read_tsv(path)
    by_block = defaultdict(list)
    by_pos = {}
    for r in rows:
        bid = r.get("scaffold_block_id", "")
        state = as_int(r.get("scaffold_phase_state", 0))
        conf = as_float(r.get("scaffold_confidence", 0.0))
        pos = as_int(r.get("pos"))
        if pos and state != 0 and conf > 0:
            by_pos[pos] = (state, conf, bid)
            if bid:
                by_block[bid].append((pos, state, conf))
    return by_block, by_pos


def scaffold_synthetic_for_region(region_positions, scaffold_by_block, scaffold_by_pos, anchor_weight, continuity_weight):
    rset = set(region_positions)
    syn = []
    # Direct scaffold sites inside region.
    for bid, items in scaffold_by_block.items():
        sub = [(p, s, c) for p, s, c in items if p in rset]
        if len(sub) >= 1:
            posmap = {p: (s, max(anchor_weight, c * 4.0)) for p, s, c in sub}
            syn.append((f"SCAFFOLD::{bid}", posmap))
    # If no direct scaffold sites are in region, create a mild continuity pseudo-fragment from nearest scaffolded internal calls
    # only when the region positions themselves were already scaffolded in by_pos.
    direct = {p: scaffold_by_pos[p] for p in region_positions if p in scaffold_by_pos}
    if direct:
        posmap = {p: (v[0], max(continuity_weight, v[1] * 3.0)) for p, v in direct.items()}
        syn.append(("SCAFFOLD::DIRECT", posmap))
    return syn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--scaffold-calls-tsv", required=True)
    ap.add_argument("--max-coverage", type=int, default=32)
    ap.add_argument("--max-reads-per-region", type=int, default=96)
    ap.add_argument("--min-obs-confidence", type=float, default=0.025)
    ap.add_argument("--min-obs-delta", type=float, default=0.45)
    ap.add_argument("--merge-logodds", type=float, default=6.0)
    ap.add_argument("--scaffold-anchor-weight", type=float, default=4.0)
    ap.add_argument("--scaffold-continuity-weight", type=float, default=2.0)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-region-summary-tsv")
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    regions = load_regions(args.regions_tsv)
    read_obs = load_obs(args.obs_tsv, args.min_obs_confidence, args.min_obs_delta)
    scaffold_blocks, scaffold_by_pos = load_scaffold_calls(args.scaffold_calls_tsv)

    rows = []; reg_summ = []
    total_fragments = 0; total_real = 0; total_syn = 0; exact_ready = 0
    frag_count_hist = Counter()

    for reg in regions:
        selected = cap_reads_for_region(reg["positions"], read_obs, args.max_coverage, args.max_reads_per_region)
        merged = merge_probabilistically(selected, args.merge_logodds)
        syn = scaffold_synthetic_for_region(reg["positions"], scaffold_blocks, scaffold_by_pos, args.scaffold_anchor_weight, args.scaffold_continuity_weight)

        frag_idx = 0; real_n = 0; syn_n = 0
        for cl in merged:
            posmap = cl["posmap"]
            # Keep singleton read fragments only in singleton regions; otherwise require bridging.
            if len(posmap) < (1 if len(reg["positions"]) == 1 else 2):
                continue
            frag_idx += 1; real_n += 1; total_real += 1
            sps = sorted(posmap)
            rows.append({
                "region_id": reg["region_id"],
                "fragment_id": f"{reg['region_id']}::R{frag_idx}",
                "source_type": "read_merged" if len(cl["members"]) > 1 else "read",
                "n_sites": len(posmap),
                "positions_json": json.dumps(sps),
                "alleles_json": json.dumps([posmap[p][0] for p in sps]),
                "weights_json": json.dumps([round(posmap[p][1], 6) for p in sps]),
                "members_json": json.dumps(cl["members"]),
                "fragment_score": round(fragment_score(posmap), 6),
            })
        for sid, posmap in syn:
            if len(posmap) < 1:
                continue
            frag_idx += 1; syn_n += 1; total_syn += 1
            sps = sorted(posmap)
            rows.append({
                "region_id": reg["region_id"],
                "fragment_id": f"{reg['region_id']}::S{frag_idx}",
                "source_type": "scaffold_synthetic",
                "n_sites": len(posmap),
                "positions_json": json.dumps(sps),
                "alleles_json": json.dumps([posmap[p][0] for p in sps]),
                "weights_json": json.dumps([round(posmap[p][1], 6) for p in sps]),
                "members_json": json.dumps([sid]),
                "fragment_score": round(fragment_score(posmap), 6),
            })
        per_site = defaultdict(int)
        for r in rows[-frag_idx:] if frag_idx else []:
            for p in json.loads(r["positions_json"]):
                per_site[p] += 1
        active_cov = max(per_site.values()) if per_site else 0
        if active_cov <= args.max_coverage and len(reg["positions"]) <= 18 and (real_n >= 1 or syn_n >= 1):
            exact_ready += 1
        total_fragments += frag_idx
        frag_count_hist[frag_idx] += 1
        reg_summ.append({
            "region_id": reg["region_id"],
            "n_sites": len(reg["positions"]),
            "real_fragments": real_n,
            "scaffold_fragments": syn_n,
            "total_fragments": frag_idx,
            "max_effective_coverage": active_cov,
        })

    fields = ["region_id","fragment_id","source_type","n_sites","positions_json","alleles_json","weights_json","members_json","fragment_score"]
    write_tsv(args.out_tsv, rows, fields)
    if args.out_region_summary_tsv:
        write_tsv(args.out_region_summary_tsv, reg_summ, ["region_id","n_sites","real_fragments","scaffold_fragments","total_fragments","max_effective_coverage"])
    summary = {
        "regions": len(regions),
        "fragments": total_fragments,
        "real_fragments": total_real,
        "scaffold_fragments": total_syn,
        "exact_ready_regions": exact_ready,
        "max_coverage": args.max_coverage,
        "merge_logodds": args.merge_logodds,
        "regions_with_no_fragments": sum(1 for r in reg_summ if int(r["total_fragments"]) == 0),
        "regions_with_real_fragments": sum(1 for r in reg_summ if int(r["real_fragments"]) > 0),
        "fragment_count_hist": dict(sorted(frag_count_hist.items())),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
