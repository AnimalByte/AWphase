#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ival(x, default=0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return int(float(x))
    except Exception:
        return default


def load_variants(path):
    obj = json.load(open(path))
    rows = obj.get("variants", obj) if isinstance(obj, dict) else obj

    out = {}
    for r in rows:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        ref = str(r.get("ref_allele", r.get("ref", ""))).upper()
        alt = str(r.get("alt_allele", r.get("alt", ""))).upper()
        if ref and alt:
            out[pos] = (ref, alt)
    return out


def load_genetic_map(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    bp = []
    cm = []

    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                cm_val = float(parts[-2])
                bp_val = int(float(parts[-1]))
            except Exception:
                continue
            bp.append(bp_val)
            cm.append(cm_val)

    if not bp:
        raise SystemExit(f"No genetic-map points loaded from {path}")

    order = np.argsort(bp)
    return (
        np.array([bp[i] for i in order], dtype=np.float64),
        np.array([cm[i] for i in order], dtype=np.float64),
    )


def interp_cm(pos, map_bp, map_cm):
    return float(np.interp(float(pos), map_bp, map_cm))


def load_phase6_blocks(local_calls_tsv):
    rows = read_tsv(local_calls_tsv)
    by_block = defaultdict(list)

    for r in rows:
        pos = ival(r.get("pos"), None)
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))

        if pos is None:
            continue
        if state == 0:
            continue
        if block in {"", "0", "unassigned"}:
            continue

        by_block[block].append((pos, state))

    for b in by_block:
        by_block[b].sort(key=lambda x: x[0])

    stats = {}
    for b, arr in by_block.items():
        poss = [p for p, _ in arr]
        stats[b] = {
            "start": min(poss),
            "end": max(poss),
            "span_bp": max(poss) - min(poss) + 1,
            "n_anchors": len(poss),
        }

    return by_block, stats


def load_truth_orientation(site_comparison_tsv, min_orientation_sites):
    rows = read_tsv(site_comparison_tsv)
    counts = defaultdict(lambda: [0, 0])  # same, opposite

    for r in rows:
        block = str(r.get("block_id", "unassigned"))
        if block in {"", "0", "unassigned"}:
            continue

        pred = ival(r.get("pred_state", r.get("local_phase_state", 0)), 0)
        truth = ival(r.get("truth_state", 0), 0)

        if pred == 0 or truth == 0:
            continue

        if pred == truth:
            counts[block][0] += 1
        elif -pred == truth:
            counts[block][1] += 1

    orient = {}
    orient_sites = {}
    orient_margin = {}

    for block, (same, opp) in counts.items():
        total = same + opp
        margin = abs(same - opp)

        if total < min_orientation_sites:
            continue
        if margin == 0:
            continue

        orient[block] = 1 if same > opp else -1
        orient_sites[block] = total
        orient_margin[block] = margin

    return orient, orient_sites, orient_margin


def load_panel_states(panel_bcf, chrom, start, end, needed_variants):
    vf = pysam.VariantFile(panel_bcf)
    samples = list(vf.header.samples)
    n_haps = len(samples) * 2

    needed = set(needed_variants)
    states = {}
    found = 0
    mismatch = 0
    nonbiallelic = 0

    for rec in vf.fetch(chrom, max(0, start - 1), end):
        pos = int(rec.pos)
        if pos not in needed:
            continue

        want_ref, want_alt = needed_variants[pos]
        ref = str(rec.ref).upper()
        alts = rec.alts or []

        if len(alts) != 1:
            nonbiallelic += 1
            continue

        alt = str(alts[0]).upper()
        if ref != want_ref or alt != want_alt:
            mismatch += 1
            continue

        alleles = np.zeros(n_haps, dtype=np.int8)
        valid = np.zeros(n_haps, dtype=bool)

        for si, sample in enumerate(samples):
            call = rec.samples[sample]
            gt = call.get("GT")
            phased = bool(call.phased)

            if gt is None or len(gt) != 2:
                continue
            if gt[0] is None or gt[1] is None:
                continue
            if not phased:
                continue
            if gt[0] not in (0, 1) or gt[1] not in (0, 1):
                continue

            alleles[2 * si] = gt[0]
            alleles[2 * si + 1] = gt[1]
            valid[2 * si] = True
            valid[2 * si + 1] = True

        states[pos] = (alleles, valid)
        found += 1

    vf.close()

    return {
        "n_haps": n_haps,
        "states": states,
        "positions_needed": len(needed),
        "positions_found": found,
        "allele_mismatch": mismatch,
        "nonbiallelic": nonbiallelic,
    }


def effective_n(weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w[w > 0]
    if w.size == 0:
        return 0.0
    return float((w.sum() ** 2) / max(float((w * w).sum()), 1e-12))


def entropy_binary(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return float(-(p * math.log(p) + (1.0 - p) * math.log(1.0 - p)))


def topk_stats(scores, k):
    vals = np.asarray(scores, dtype=np.float64)
    vals = vals[vals > 0]
    if vals.size == 0:
        return 0.0, 0.0, 0.0, 0.0

    vals = np.sort(vals)[::-1]
    top = vals[:k]
    top1 = float(top[0]) if top.size > 0 else 0.0
    top2 = float(top[1]) if top.size > 1 else 0.0
    return float(top.sum()), top1, top2, effective_n(top)


def block_scores(edge_anchors, panel, map_bp, map_cm, midpoint_cm, decay_cm):
    n_haps = panel["n_haps"]
    plus = np.zeros(n_haps, dtype=np.float64)
    minus = np.zeros(n_haps, dtype=np.float64)

    usable = 0
    positions_found = 0

    for pos, state in edge_anchors:
        if pos not in panel["states"]:
            continue

        alleles, valid = panel["states"][pos]
        cm = interp_cm(pos, map_bp, map_cm)
        dcm = abs(cm - midpoint_cm)
        w = math.exp(-dcm / max(decay_cm, 1e-9))

        plus_expected = 1 if state == 1 else 0
        minus_expected = 1 - plus_expected

        plus_match = valid & (alleles == plus_expected)
        minus_match = valid & (alleles == minus_expected)

        plus[plus_match] += w
        minus[minus_match] += w

        usable += 1
        positions_found += 1

    return plus, minus, usable, positions_found


def score_bridge(a_edge, b_edge, panel, map_bp, map_cm, decay_cm, top_k):
    all_pos = [p for p, _ in a_edge] + [p for p, _ in b_edge]
    mid_bp = (min(all_pos) + max(all_pos)) / 2.0
    mid_cm = interp_cm(mid_bp, map_bp, map_cm)

    a_plus, a_minus, a_usable, _ = block_scores(a_edge, panel, map_bp, map_cm, mid_cm, decay_cm)
    b_plus, b_minus, b_usable, _ = block_scores(b_edge, panel, map_bp, map_cm, mid_cm, decay_cm)

    same_per_hap = np.minimum(a_plus, b_plus) + np.minimum(a_minus, b_minus)
    flip_per_hap = np.minimum(a_plus, b_minus) + np.minimum(a_minus, b_plus)

    same_support, same_top1, same_top2, same_effn = topk_stats(same_per_hap, top_k)
    flip_support, flip_top1, flip_top2, flip_effn = topk_stats(flip_per_hap, top_k)

    total = same_support + flip_support

    if total > 0:
        pred_relation = 1 if same_support >= flip_support else -1
        confidence = max(same_support, flip_support) / total
        margin = abs(same_support - flip_support)
        entropy = entropy_binary(confidence)
        log_margin = math.log((same_support + 1.0) / (flip_support + 1.0))
    else:
        pred_relation = 0
        confidence = 0.0
        margin = 0.0
        entropy = 0.0
        log_margin = 0.0

    return {
        "bridge_pred_relation": pred_relation,
        "bridge_a_usable_anchors": a_usable,
        "bridge_b_usable_anchors": b_usable,
        "bridge_same_support": same_support,
        "bridge_flip_support": flip_support,
        "bridge_total_support": total,
        "bridge_margin": margin,
        "bridge_confidence": confidence,
        "bridge_entropy": entropy,
        "bridge_log_margin_same_over_flip": log_margin,
        "bridge_abs_log_margin": abs(log_margin),
        "bridge_same_top1": same_top1,
        "bridge_same_top2": same_top2,
        "bridge_flip_top1": flip_top1,
        "bridge_flip_top2": flip_top2,
        "bridge_top1_max": max(same_top1, flip_top1),
        "bridge_top2_max": max(same_top2, flip_top2),
        "bridge_top_margin": max(same_top1, flip_top1) - max(same_top2, flip_top2),
        "bridge_effective_donors": max(same_effn, flip_effn),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--truth-comparison-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--genetic-map", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--max-gap-bp", type=int, default=100000)
    ap.add_argument("--max-next-blocks", type=int, default=3)
    ap.add_argument("--max-anchors-each-block", type=int, default=12)
    ap.add_argument("--min-anchors-each-block", type=int, default=2)
    ap.add_argument("--min-orientation-sites", type=int, default=2)
    ap.add_argument("--decay-cm", type=float, default=0.05)
    ap.add_argument("--top-k", type=int, default=64)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    variants = load_variants(args.variant_json)
    map_bp, map_cm = load_genetic_map(args.genetic_map)

    by_block, stats = load_phase6_blocks(args.local_calls_tsv)
    orient, orient_sites, orient_margin = load_truth_orientation(
        args.truth_comparison_tsv,
        args.min_orientation_sites,
    )

    blocks = [
        b for b in by_block
        if b in stats and len(by_block[b]) >= args.min_anchors_each_block
    ]
    blocks.sort(key=lambda b: (stats[b]["start"], stats[b]["end"]))

    raw_candidates = []
    needed_positions = set()

    for i, a in enumerate(blocks):
        for j in range(i + 1, min(len(blocks), i + 1 + args.max_next_blocks)):
            b = blocks[j]
            gap = stats[b]["start"] - stats[a]["end"]
            if gap < 0:
                continue
            if gap > args.max_gap_bp:
                break

            a_edge = by_block[a][-args.max_anchors_each_block:]
            b_edge = by_block[b][:args.max_anchors_each_block]

            if len(a_edge) < args.min_anchors_each_block:
                continue
            if len(b_edge) < args.min_anchors_each_block:
                continue

            for pos, _ in a_edge + b_edge:
                if pos in variants:
                    needed_positions.add(pos)

            raw_candidates.append((a, b, gap, a_edge, b_edge))

    needed_variants = {p: variants[p] for p in needed_positions if p in variants}

    panel = load_panel_states(
        args.panel_bcf,
        args.chrom,
        args.start,
        args.end,
        needed_variants,
    )

    rows = []

    for a, b, gap, a_edge, b_edge in raw_candidates:
        if a not in orient or b not in orient:
            continue

        true_relation = 1 if orient[a] == orient[b] else -1

        score = score_bridge(
            a_edge=a_edge,
            b_edge=b_edge,
            panel=panel,
            map_bp=map_bp,
            map_cm=map_cm,
            decay_cm=args.decay_cm,
            top_k=args.top_k,
        )

        pred_relation = int(score["bridge_pred_relation"])
        if pred_relation == 0:
            continue

        label = 1 if pred_relation == true_relation else 0

        a_end_cm = interp_cm(stats[a]["end"], map_bp, map_cm)
        b_start_cm = interp_cm(stats[b]["start"], map_bp, map_cm)
        gap_cm = max(0.0, b_start_cm - a_end_cm)

        row = {
            "source_window": args.source_window,
            "label": label,
            "block_a": a,
            "block_b": b,
            "true_relation": true_relation,
            "pred_relation": pred_relation,
            "gap_bp": gap,
            "gap_cm": gap_cm,
            "block_a_start": stats[a]["start"],
            "block_a_end": stats[a]["end"],
            "block_a_span_bp": stats[a]["span_bp"],
            "block_a_n_anchors": stats[a]["n_anchors"],
            "block_b_start": stats[b]["start"],
            "block_b_end": stats[b]["end"],
            "block_b_span_bp": stats[b]["span_bp"],
            "block_b_n_anchors": stats[b]["n_anchors"],
            "bridge_anchor_count_a": len(a_edge),
            "bridge_anchor_count_b": len(b_edge),
            "block_a_orientation_sites": orient_sites.get(a, 0),
            "block_b_orientation_sites": orient_sites.get(b, 0),
            "block_a_orientation_margin": orient_margin.get(a, 0),
            "block_b_orientation_margin": orient_margin.get(b, 0),
        }
        row.update(score)
        rows.append(row)

    fields = []
    seen = set()
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    pos = sum(int(r["label"]) for r in rows)
    neg = len(rows) - pos

    print(json.dumps({
        "source_window": args.source_window,
        "blocks": len(blocks),
        "raw_candidates": len(raw_candidates),
        "labeled_rows": len(rows),
        "positives": pos,
        "negatives": neg,
        "panel_positions_needed": panel["positions_needed"],
        "panel_positions_found": panel["positions_found"],
        "panel_allele_mismatch": panel["allele_mismatch"],
        "panel_nonbiallelic": panel["nonbiallelic"],
        "out": args.out_tsv,
    }, indent=2))


if __name__ == "__main__":
    main()
