#!/usr/bin/env python3
import argparse
import csv
import math
import sys
from pathlib import Path
from bisect import bisect_left
from collections import defaultdict

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def fval(x, default=0.0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return float(x)
    except Exception:
        return default

def ival(x, default=0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return int(float(x))
    except Exception:
        return default

def safe_log1p(x):
    x = fval(x, 0.0)
    if x < 0:
        x = 0.0
    return math.log1p(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    train = read_tsv(args.training_tsv)
    calls = read_tsv(args.local_calls_phase6c_tsv)

    anchors = []
    anchors_by_block = defaultdict(list)

    for r in calls:
        pos = ival(r.get("pos"), None)
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))

        if pos is None:
            continue

        if state != 0 and block not in {"", "0", "unassigned"}:
            anchors.append(pos)
            anchors_by_block[block].append(pos)

    anchors.sort()
    for b in anchors_by_block:
        anchors_by_block[b].sort()

    block_stats = {}
    for block, poss in anchors_by_block.items():
        if poss:
            span = max(poss) - min(poss) + 1
            block_stats[block] = {
                "phase6_block_anchor_count": len(poss),
                "phase6_block_span_bp": span,
                "phase6_block_density_per_kb": len(poss) / max(span / 1000.0, 1e-9),
            }

    def nearest_dist(pos, arr):
        if not arr:
            return -1
        i = bisect_left(arr, pos)
        dists = []
        if i < len(arr):
            dists.append(abs(arr[i] - pos))
        if i > 0:
            dists.append(abs(arr[i - 1] - pos))
        return min(dists) if dists else -1

    def side_features(pos, arr, radius):
        if not arr:
            return (0, 0, -1, -1)

        i = bisect_left(arr, pos)

        left_count = 0
        j = i - 1
        left_nearest = -1
        while j >= 0 and pos - arr[j] <= radius:
            if left_nearest < 0:
                left_nearest = pos - arr[j]
            left_count += 1
            j -= 1

        right_count = 0
        j = i
        right_nearest = -1
        while j < len(arr) and arr[j] - pos <= radius:
            if right_nearest < 0:
                right_nearest = arr[j] - pos
            right_count += 1
            j += 1

        return left_count, right_count, left_nearest, right_nearest

    out = []

    for r in train:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        block = str(rr.get("block_id", "unassigned"))

        panel_confidence = fval(rr.get("panel_confidence"))
        panel_margin = fval(rr.get("panel_margin"))
        panel_support = fval(rr.get("panel_support"))
        panel_conflict = fval(rr.get("panel_conflict"))
        panel_haplotypes = fval(rr.get("panel_haplotypes"))
        anchors_n = fval(rr.get("anchors"))
        best_second = fval(rr.get("best_vs_second_margin"))

        rr["log_panel_margin"] = safe_log1p(panel_margin)
        rr["log_panel_support"] = safe_log1p(panel_support)
        rr["log_panel_conflict"] = safe_log1p(panel_conflict)
        rr["log_best_vs_second_margin"] = safe_log1p(best_second)

        rr["panel_margin_per_haplotype"] = panel_margin / max(panel_haplotypes, 1.0)
        rr["panel_support_per_haplotype"] = panel_support / max(panel_haplotypes, 1.0)
        rr["panel_conflict_per_haplotype"] = panel_conflict / max(panel_haplotypes, 1.0)
        rr["panel_haplotypes_per_anchor"] = panel_haplotypes / max(anchors_n, 1.0)
        rr["confidence_x_log_margin"] = panel_confidence * safe_log1p(panel_margin)
        rr["confidence_x_log_support_conflict_ratio"] = panel_confidence * safe_log1p(rr.get("support_conflict_ratio"))

        if pos is None:
            rr["nearest_phase6_anchor_dist_bp"] = -1
            rr["nearest_same_block_anchor_dist_bp"] = -1
            rr["left_anchors_20kb"] = 0
            rr["right_anchors_20kb"] = 0
            rr["two_sided_anchors_20kb"] = 0
            rr["left_nearest_anchor_dist_bp"] = -1
            rr["right_nearest_anchor_dist_bp"] = -1
        else:
            same_block = anchors_by_block.get(block, [])

            rr["nearest_phase6_anchor_dist_bp"] = nearest_dist(pos, anchors)
            rr["nearest_same_block_anchor_dist_bp"] = nearest_dist(pos, same_block)

            l20, r20, ld20, rd20 = side_features(pos, same_block, 20000)
            l50, r50, ld50, rd50 = side_features(pos, same_block, 50000)

            rr["left_anchors_20kb"] = l20
            rr["right_anchors_20kb"] = r20
            rr["two_sided_anchors_20kb"] = 1 if l20 > 0 and r20 > 0 else 0
            rr["left_nearest_anchor_dist_bp"] = ld20
            rr["right_nearest_anchor_dist_bp"] = rd20

            rr["left_anchors_50kb"] = l50
            rr["right_anchors_50kb"] = r50
            rr["two_sided_anchors_50kb"] = 1 if l50 > 0 and r50 > 0 else 0

        bs = block_stats.get(block, {})
        rr["phase6_block_anchor_count"] = bs.get("phase6_block_anchor_count", 0)
        rr["phase6_block_span_bp"] = bs.get("phase6_block_span_bp", 0)
        rr["phase6_block_density_per_kb"] = bs.get("phase6_block_density_per_kb", 0.0)
        rr["log_phase6_block_span_bp"] = safe_log1p(rr["phase6_block_span_bp"])
        rr["log_phase6_block_anchor_count"] = safe_log1p(rr["phase6_block_anchor_count"])

        out.append(rr)

    fields = []
    seen = set()
    for r in out:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print({
        "input_rows": len(train),
        "output_rows": len(out),
        "anchors": len(anchors),
        "blocks": len(anchors_by_block),
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
