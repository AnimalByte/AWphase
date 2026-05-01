#!/usr/bin/env python3
import argparse
import csv
import json
import bisect
from pathlib import Path
from collections import defaultdict
import pysam

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            return ival(row.get(k), 0)
    return 0

def set_state(row, s):
    if "local_phase_state" in row:
        row["local_phase_state"] = int(s)
    if "phase_state" in row:
        row["phase_state"] = int(s)

def block_id(row):
    b = str(row.get("block_id", "")).strip()
    if b == "" or b.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return b

def load_variants(path, chrom, start, end):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    variants = {}
    order = []

    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue

        c = str(r.get("chrom", r.get("contig", chrom or "")))
        if c and c != chrom:
            continue
        if pos < start or pos > end:
            continue

        ref = str(r.get("ref_allele", r.get("ref", ""))).upper()
        alt = str(r.get("alt_allele", r.get("alt", ""))).upper()

        if len(ref) != 1 or len(alt) != 1:
            continue

        variants[pos] = {"pos": pos, "ref": ref, "alt": alt}
        order.append(pos)

    return variants, sorted(set(order))

def sample_hap_alleles(sample_data):
    gt = sample_data.get("GT")
    phased = getattr(sample_data, "phased", False)

    if not phased or gt is None or len(gt) != 2 or None in gt:
        return None

    a, b = gt
    if a not in (0, 1) or b not in (0, 1):
        return None

    return int(a), int(b)

def load_panel_haplotypes_stream(panel_bcf, chrom, variants, needed_positions, start, end):
    needed = set(needed_positions)
    states_by_pos = {}
    found = set()
    allele_mismatch = 0
    nonbiallelic = 0
    retained_states = 0

    vf = pysam.VariantFile(panel_bcf)

    for rec in vf.fetch(chrom, start - 1, end):
        pos = int(rec.pos)
        if pos not in needed:
            continue

        v = variants.get(pos)
        if v is None:
            continue

        if rec.alts is None or len(rec.alts) != 1:
            nonbiallelic += 1
            continue

        if str(rec.ref).upper() != v["ref"] or str(rec.alts[0]).upper() != v["alt"]:
            allele_mismatch += 1
            continue

        d = {}
        for sample in rec.samples:
            h = sample_hap_alleles(rec.samples[sample])
            if h is not None:
                d[sample] = h

        states_by_pos[pos] = d
        found.add(pos)
        retained_states += len(d)

    vf.close()

    for pos in needed:
        states_by_pos.setdefault(pos, {})

    return states_by_pos, {
        "needed_panel_positions": len(needed),
        "panel_positions_found": len(found),
        "panel_positions_missing": len(needed - found),
        "panel_allele_mismatch_records": allele_mismatch,
        "panel_nonbiallelic_records": nonbiallelic,
        "retained_sample_states": retained_states,
    }

def nearby_anchors(pos, anchor_positions, max_dist, max_anchors_each_side):
    i = bisect.bisect_left(anchor_positions, pos)
    out = []

    j = i - 1
    n = 0
    while j >= 0 and n < max_anchors_each_side:
        d = pos - anchor_positions[j]
        if d > max_dist:
            break
        out.append((anchor_positions[j], d))
        j -= 1
        n += 1

    j = i
    n = 0
    while j < len(anchor_positions) and n < max_anchors_each_side:
        d = anchor_positions[j] - pos
        if d > max_dist:
            break
        out.append((anchor_positions[j], d))
        j += 1
        n += 1

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--backbone-local-calls-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fill-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--max-anchor-dist-bp", type=int, default=5000)
    ap.add_argument("--max-anchors-each-side", type=int, default=4)
    ap.add_argument("--min-panel-samples", type=int, default=25)
    ap.add_argument("--min-confidence", type=float, default=0.90)
    ap.add_argument("--min-margin", type=float, default=10.0)
    ap.add_argument("--min-best-vs-second-margin", type=float, default=5.0)
    args = ap.parse_args()

    variants, variant_order = load_variants(args.variant_json, args.chrom, args.start, args.end)
    rows = read_tsv(args.backbone_local_calls_tsv)

    by_pos = {}
    anchors = {}
    anchor_positions = []

    for r in rows:
        pos = ival(r.get("pos"), 0)
        if pos <= 0:
            continue
        by_pos[pos] = r

        s = state(r)
        b = block_id(r)
        if s != 0 and b and pos in variants:
            anchors[pos] = {"state": s, "block_id": b}
            anchor_positions.append(pos)

    anchor_positions = sorted(set(anchor_positions))

    candidate_positions = [
        p for p in variant_order
        if p in by_pos and state(by_pos[p]) == 0
    ]

    needed = set(candidate_positions)
    for p in candidate_positions:
        for a, _ in nearby_anchors(p, anchor_positions, args.max_anchor_dist_bp, args.max_anchors_each_side):
            needed.add(a)

    panel_haps, panel_summary = load_panel_haplotypes_stream(
        args.panel_bcf, args.chrom, variants, sorted(needed), args.start, args.end
    )

    fill_rows = []
    chosen = {}

    candidates_with_anchor = 0
    candidates_with_panel = 0
    candidates_with_votes = 0

    for pos in candidate_positions:
        target_haps = panel_haps.get(pos, {})
        if target_haps:
            candidates_with_panel += 1

        anchors_here = nearby_anchors(pos, anchor_positions, args.max_anchor_dist_bp, args.max_anchors_each_side)
        if anchors_here:
            candidates_with_anchor += 1

        block_votes = defaultdict(lambda: {
            "plus": 0.0,
            "minus": 0.0,
            "panel_haplotypes": 0,
            "samples": set(),
            "anchors": set(),
        })

        for anchor_pos, dist in anchors_here:
            anchor_info = anchors[anchor_pos]
            anchor_haps = panel_haps.get(anchor_pos, {})
            if not anchor_haps or not target_haps:
                continue

            shared = set(anchor_haps) & set(target_haps)
            if len(shared) < args.min_panel_samples:
                continue

            dist_weight = 1.0 / (1.0 + dist / 5000.0)

            for sample in shared:
                ah = anchor_haps[sample]
                th = target_haps[sample]

                # Use both haplotypes from each phased panel sample.
                for hap_i in (0, 1):
                    anchor_allele = ah[hap_i]
                    target_allele = th[hap_i]

                    # Only informative if target and anchor are both variant-discriminating.
                    # allele state: REF=-1, ALT=+1
                    anchor_state = 1 if anchor_allele == 1 else -1
                    target_state = 1 if target_allele == 1 else -1

                    relation = target_state * anchor_state
                    pred = anchor_info["state"] * relation
                    block = anchor_info["block_id"]

                    if pred == 1:
                        block_votes[block]["plus"] += dist_weight
                    else:
                        block_votes[block]["minus"] += dist_weight

                    block_votes[block]["panel_haplotypes"] += 1
                    block_votes[block]["samples"].add(sample)
                    block_votes[block]["anchors"].add(anchor_pos)

        block_candidates = []

        for block, v in block_votes.items():
            plus = v["plus"]
            minus = v["minus"]
            support = max(plus, minus)
            conflict = min(plus, minus)
            margin = support - conflict
            total = support + conflict
            conf = support / total if total else 0.0
            pred = 1 if plus >= minus else -1

            block_candidates.append({
                "pos": pos,
                "block_id": block,
                "pred_state": pred,
                "support": support,
                "conflict": conflict,
                "margin": margin,
                "confidence": conf,
                "panel_samples": len(v["samples"]),
                "panel_haplotypes": v["panel_haplotypes"],
                "anchors": len(v["anchors"]),
            })

        if block_candidates:
            candidates_with_votes += 1
        else:
            continue

        block_candidates.sort(
            key=lambda x: (x["margin"], x["confidence"], x["panel_haplotypes"], x["anchors"]),
            reverse=True,
        )

        best = block_candidates[0]
        second_margin = block_candidates[1]["margin"] if len(block_candidates) > 1 else 0.0
        best_vs_second = best["margin"] - second_margin

        accepted = (
            best["panel_samples"] >= args.min_panel_samples
            and best["confidence"] >= args.min_confidence
            and best["margin"] >= args.min_margin
            and best_vs_second >= args.min_best_vs_second_margin
        )

        best["best_vs_second_margin"] = best_vs_second
        best["n_block_candidates"] = len(block_candidates)
        best["accepted"] = int(accepted)

        fill_rows.append(best)

        if accepted:
            chosen[pos] = best

    out_rows = []
    filled = 0

    for r in rows:
        rr = dict(r)
        pos = ival(rr.get("pos"), 0)

        if pos in chosen and state(rr) == 0:
            f = chosen[pos]
            set_state(rr, f["pred_state"])
            rr["block_id"] = f["block_id"]
            rr["phase7a_panel_filled"] = 1
            rr["phase7a_panel_confidence"] = f"{f['confidence']:.6f}"
            rr["phase7a_panel_margin"] = f"{f['margin']:.6f}"
            rr["phase7a_panel_samples"] = f["panel_samples"]
            rr["phase7a_panel_haplotypes"] = f["panel_haplotypes"]
            rr["phase7a_panel_anchors"] = f["anchors"]
            filled += 1
        else:
            rr.setdefault("phase7a_panel_filled", 0)
            rr.setdefault("phase7a_panel_confidence", "")
            rr.setdefault("phase7a_panel_margin", "")
            rr.setdefault("phase7a_panel_samples", "")
            rr.setdefault("phase7a_panel_haplotypes", "")
            rr.setdefault("phase7a_panel_anchors", "")

        out_rows.append(rr)

    write_tsv(args.out_tsv, out_rows)

    fill_fields = [
        "pos", "block_id", "pred_state", "accepted",
        "support", "conflict", "margin", "confidence",
        "panel_samples", "panel_haplotypes", "anchors",
        "best_vs_second_margin", "n_block_candidates",
    ]

    Path(args.out_fill_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_fill_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fill_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(fill_rows)

    summary = {
        "panel_bcf": args.panel_bcf,
        "backbone_local_calls_tsv": args.backbone_local_calls_tsv,
        "variant_json": args.variant_json,
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "backbone_anchors": len(anchor_positions),
        "candidate_unphased_positions": len(candidate_positions),
        "candidates_with_nearby_anchor": candidates_with_anchor,
        "candidates_with_panel_state": candidates_with_panel,
        "candidates_with_panel_votes": candidates_with_votes,
        "fill_candidate_sites": len(fill_rows),
        "filled_sites": filled,
        "max_anchor_dist_bp": args.max_anchor_dist_bp,
        "max_anchors_each_side": args.max_anchors_each_side,
        "min_panel_samples": args.min_panel_samples,
        "min_confidence": args.min_confidence,
        "min_margin": args.min_margin,
        "min_best_vs_second_margin": args.min_best_vs_second_margin,
        **panel_summary,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
