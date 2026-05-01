#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import defaultdict, Counter

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
    return 0

def block_id(row):
    b = str(row.get("block_id", "")).strip()
    if b == "" or b.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return b

def load_calls(path):
    rows = {}
    blocks = defaultdict(list)

    for r in read_tsv(path):
        try:
            pos = int(float(r["pos"]))
        except Exception:
            continue

        s = state(r)
        b = block_id(r)

        rows[pos] = {
            "row": r,
            "state": s,
            "block_id": b,
        }

        if s != 0 and b:
            blocks[b].append(pos)

    for b in blocks:
        blocks[b].sort()

    return rows, blocks

def load_recall_gap(path):
    rows = []
    for r in read_tsv(path):
        try:
            r["pos"] = int(float(r["pos"]))
            r["reads_covering"] = int(float(r.get("reads_covering", 0)))
            r["templates_covering"] = int(float(r.get("templates_covering", 0)))
        except Exception:
            continue
        rows.append(r)
    return rows

def nearest_distance(pos, sorted_positions):
    if not sorted_positions:
        return None

    import bisect
    i = bisect.bisect_left(sorted_positions, pos)
    best = None

    for j in (i - 1, i):
        if 0 <= j < len(sorted_positions):
            d = abs(pos - sorted_positions[j])
            best = d if best is None else min(best, d)

    return best

def bucket_block_size(n):
    if n <= 0:
        return "0"
    if n == 1:
        return "1_singleton"
    if n == 2:
        return "2"
    if n <= 5:
        return "3_5"
    if n <= 10:
        return "6_10"
    if n <= 25:
        return "11_25"
    if n <= 100:
        return "26_100"
    return "gt100"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6c-tsv", required=True)
    ap.add_argument("--whatshap-tsv", required=True)
    ap.add_argument("--whatshap-only-sites-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    p6, p6_blocks = load_calls(args.phase6c_tsv)
    wh, wh_blocks = load_calls(args.whatshap_tsv)
    gap = load_recall_gap(args.whatshap_only_sites_tsv)

    p6_nonzero_positions = sorted([p for p, x in p6.items() if x["state"] != 0])
    wh_nonzero_positions = sorted([p for p, x in wh.items() if x["state"] != 0])

    out_rows = []

    bucket_counts = Counter()
    type_counts = Counter()
    near_p6_counts = Counter()

    for r in gap:
        pos = r["pos"]
        wh_info = wh.get(pos, {})
        b = wh_info.get("block_id", "")
        block_positions = wh_blocks.get(b, [])
        block_n = len(block_positions)
        block_span = block_positions[-1] - block_positions[0] + 1 if block_positions else 0

        idx = block_positions.index(pos) if pos in block_positions else -1
        left_dist = pos - block_positions[idx - 1] if idx > 0 else ""
        right_dist = block_positions[idx + 1] - pos if idx >= 0 and idx + 1 < len(block_positions) else ""

        nearest_p6 = nearest_distance(pos, p6_nonzero_positions)
        nearest_wh = nearest_distance(pos, wh_nonzero_positions)

        if nearest_p6 is None:
            near_bucket = "none"
        elif nearest_p6 <= 500:
            near_bucket = "le500"
        elif nearest_p6 <= 1000:
            near_bucket = "le1000"
        elif nearest_p6 <= 5000:
            near_bucket = "le5000"
        else:
            near_bucket = "gt5000"

        bucket = bucket_block_size(block_n)

        bucket_counts[bucket] += 1
        type_counts[r.get("type", "unknown")] += 1
        near_p6_counts[near_bucket] += 1

        rr = dict(r)
        rr.update({
            "whatshap_block_id": b,
            "whatshap_block_n_sites": block_n,
            "whatshap_block_span_bp": block_span,
            "left_wh_phased_neighbor_dist": left_dist,
            "right_wh_phased_neighbor_dist": right_dist,
            "nearest_phase6c_phased_dist": "" if nearest_p6 is None else nearest_p6,
            "nearest_whatshap_phased_dist": "" if nearest_wh is None else nearest_wh,
            "whatshap_block_size_bucket": bucket,
            "nearest_phase6c_bucket": near_bucket,
        })

        out_rows.append(rr)

    fields = [
        "pos", "ref", "alt", "type",
        "whatshap_state", "phase6c_state",
        "reads_covering", "templates_covering",
        "whatshap_block_id", "whatshap_block_n_sites", "whatshap_block_span_bp",
        "left_wh_phased_neighbor_dist", "right_wh_phased_neighbor_dist",
        "nearest_phase6c_phased_dist", "nearest_whatshap_phased_dist",
        "whatshap_block_size_bucket", "nearest_phase6c_bucket",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    summary = {
        "whatshap_only_sites": len(out_rows),
        "type_counts": dict(type_counts),
        "whatshap_block_size_bucket_counts": dict(bucket_counts),
        "nearest_phase6c_bucket_counts": dict(near_p6_counts),
        "whatshap_only_in_singleton_blocks": bucket_counts.get("1_singleton", 0),
        "whatshap_only_in_blocks_gt100": bucket_counts.get("gt100", 0),
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
