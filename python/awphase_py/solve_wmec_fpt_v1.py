#!/usr/bin/env python3
import argparse, csv, itertools, json
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def mismatch_cost(assign, frag):
    # orientation-free fragment: it can come from hap1 or hap2
    c1 = 0.0; c2 = 0.0
    for pos, allele, weight in frag:
        s = assign[pos]
        c1 += 0.0 if allele == s else weight
        c2 += 0.0 if allele == -s else weight
    return min(c1, c2)


def prior_penalty(assign, positions, site_meta, left_anchor_state, right_anchor_state,
                  donor_lambda, ranker_lambda, anchor_lambda, switch_lambda):
    pen = 0.0
    for m in site_meta:
        p = int(m["pos"]); s = assign[p]
        donor = float(m.get("donor_bias", 0.0))
        top_state = int(float(m.get("top_state", 0)))
        margin = float(m.get("margin", 0.0))
        if donor != 0.0 and s != (1 if donor > 0 else -1):
            pen += donor_lambda * min(abs(donor), 1.0)
        if top_state in (-1, 1) and s != top_state:
            pen += ranker_lambda * min(max(margin, 0.0), 2.0)
    if left_anchor_state in (-1, 1) and assign[positions[0]] != left_anchor_state:
        pen += anchor_lambda
    if right_anchor_state in (-1, 1) and assign[positions[-1]] != right_anchor_state:
        pen += anchor_lambda
    switches = sum(1 for a, b in zip(positions, positions[1:]) if assign[a] != assign[b])
    pen += switch_lambda * switches
    return pen


def total_cost(assign, positions, frags, site_meta, left_anchor_state, right_anchor_state,
               donor_lambda, ranker_lambda, anchor_lambda, switch_lambda):
    read_cost = sum(mismatch_cost(assign, frag) for frag in frags)
    pri = prior_penalty(assign, positions, site_meta, left_anchor_state, right_anchor_state,
                        donor_lambda, ranker_lambda, anchor_lambda, switch_lambda)
    return read_cost + pri, read_cost, pri


def exact_solve(positions, frags, site_meta, left_anchor_state, right_anchor_state,
                donor_lambda, ranker_lambda, anchor_lambda, switch_lambda):
    best = None; second = None
    for bits in itertools.product((1, -1), repeat=len(positions)):
        assign = {p: s for p, s in zip(positions, bits)}
        cost, read_cost, pri = total_cost(assign, positions, frags, site_meta, left_anchor_state, right_anchor_state,
                                          donor_lambda, ranker_lambda, anchor_lambda, switch_lambda)
        rec = (cost, read_cost, pri, assign)
        if best is None or cost < best[0]:
            second = best; best = rec
        elif second is None or cost < second[0]:
            second = rec
    margin = (second[0] - best[0]) if second is not None else 999.0
    return best[3], best[0], best[1], best[2], margin, True


def beam_solve(positions, frags, site_meta, left_anchor_state, right_anchor_state,
               donor_lambda, ranker_lambda, anchor_lambda, switch_lambda, beam_width):
    beam = [(0.0, {})]
    for i, p in enumerate(positions):
        new = []
        for cost, partial in beam:
            for s in (1, -1):
                nxt = dict(partial); nxt[p] = s
                approx = cost
                m = site_meta[i]
                donor = float(m.get("donor_bias", 0.0))
                top_state = int(float(m.get("top_state", 0)))
                margin = float(m.get("margin", 0.0))
                if donor != 0.0 and s != (1 if donor > 0 else -1):
                    approx += donor_lambda * min(abs(donor), 1.0)
                if top_state in (-1, 1) and s != top_state:
                    approx += ranker_lambda * min(max(margin, 0.0), 2.0)
                if i == 0 and left_anchor_state in (-1,1) and s != left_anchor_state:
                    approx += anchor_lambda
                if i > 0 and nxt[positions[i-1]] != s:
                    approx += switch_lambda
                new.append((approx, nxt))
        new.sort(key=lambda x: x[0])
        beam = new[:beam_width]
    rescored = []
    for _, assign in beam:
        cost, read_cost, pri = total_cost(assign, positions, frags, site_meta, left_anchor_state, right_anchor_state,
                                          donor_lambda, ranker_lambda, anchor_lambda, switch_lambda)
        rescored.append((cost, read_cost, pri, assign))
    rescored.sort(key=lambda x: x[0])
    best = rescored[0]; second = rescored[1] if len(rescored) > 1 else None
    margin = (second[0] - best[0]) if second is not None else 999.0
    return best[3], best[0], best[1], best[2], margin, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wmec-matrix-tsv", required=True)
    ap.add_argument("--max-exact-sites", type=int, default=18)
    ap.add_argument("--beam-width", type=int, default=256)
    ap.add_argument("--min-real-fragments", type=int, default=1)
    ap.add_argument("--allow-scaffold-only", action='store_true')
    ap.add_argument("--donor-lambda", type=float, default=0.05)
    ap.add_argument("--ranker-lambda", type=float, default=0.04)
    ap.add_argument("--anchor-lambda", type=float, default=0.75)
    ap.add_argument("--switch-lambda", type=float, default=0.10)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-region-diagnostics-tsv")
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    rows_in = read_tsv(args.wmec_matrix_tsv)
    out = []; diag = []
    exact_regions = 0; total_regions = 0; skipped = Counter()
    for row in rows_in:
        total_regions += 1
        rid = row["region_id"]
        positions = [int(x) for x in json.loads(row["positions_json"])]
        site_meta = json.loads(row["site_meta_json"])
        fr_raw = json.loads(row["fragments_json"])
        frags = []
        real_fragments = 0; scaffold_fragments = 0
        for f in fr_raw:
            ps = json.loads(f["positions_json"])
            alleles = json.loads(f["alleles_json"])
            weights = json.loads(f["weights_json"])
            frag = [(int(p), int(a), float(w)) for p, a, w in zip(ps, alleles, weights) if int(p) in set(positions)]
            if not frag: continue
            frags.append(frag)
            if str(f.get("source_type", "")).startswith("scaffold"):
                scaffold_fragments += 1
            else:
                real_fragments += 1
        reason = "solved"
        if not positions:
            skipped["no_positions"] += 1; continue
        if not frags:
            skipped["no_fragments"] += 1; continue
        if real_fragments < args.min_real_fragments and not (args.allow_scaffold_only and scaffold_fragments > 0):
            skipped["too_few_real_fragments"] += 1; continue
        left = int(float(row.get("left_anchor_state", 0)))
        right = int(float(row.get("right_anchor_state", 0)))
        if len(positions) <= args.max_exact_sites:
            assign, obj, read_obj, pri_obj, margin, used_exact = exact_solve(
                positions, frags, site_meta, left, right,
                args.donor_lambda, args.ranker_lambda, args.anchor_lambda, args.switch_lambda)
            exact_regions += 1
        else:
            assign, obj, read_obj, pri_obj, margin, used_exact = beam_solve(
                positions, frags, site_meta, left, right,
                args.donor_lambda, args.ranker_lambda, args.anchor_lambda, args.switch_lambda,
                args.beam_width)
        avg_donor = sum(float(m.get("donor_bias", 0.0)) for m in site_meta)/max(1, len(site_meta))
        anchor_ok = int((left in (0, assign[positions[0]])) and (right in (0, assign[positions[-1]])))
        switches = sum(1 for a,b in zip(positions, positions[1:]) if assign[a] != assign[b])
        diag.append({
            "region_id": rid,
            "n_sites": len(positions),
            "n_fragments": len(frags),
            "real_fragments": real_fragments,
            "scaffold_fragments": scaffold_fragments,
            "used_exact": int(used_exact),
            "objective": round(obj, 6),
            "read_objective": round(read_obj, 6),
            "prior_objective": round(pri_obj, 6),
            "margin": round(margin, 6),
            "anchor_ok": anchor_ok,
            "switches": switches,
            "avg_donor_bias": round(avg_donor, 6),
            "reason": reason,
        })
        for p in positions:
            out.append({
                "region_id": rid,
                "pos": p,
                "resolved_state": assign[p],
                "region_objective_best": round(obj, 6),
                "region_margin": round(margin, 6),
                "used_exact": int(used_exact),
                "anchor_ok": anchor_ok,
                "avg_donor_bias": round(avg_donor, 6),
                "real_fragments": real_fragments,
                "scaffold_fragments": scaffold_fragments,
            })

    write_tsv(args.out_tsv, out, list(out[0].keys()) if out else ["region_id","pos","resolved_state"])
    if args.out_region_diagnostics_tsv:
        write_tsv(args.out_region_diagnostics_tsv, diag, list(diag[0].keys()) if diag else ["region_id"])
    summary = {
        "regions": total_regions,
        "resolved_positions": len(out),
        "exact_regions": exact_regions,
        "max_exact_sites": args.max_exact_sites,
        "solved_regions": len(diag),
        "skipped": dict(skipped),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
