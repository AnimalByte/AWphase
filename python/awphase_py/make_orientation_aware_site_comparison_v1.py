#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)

def ival(x, default=0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return int(float(x))
    except Exception:
        return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-comparison-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-tsv", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.site_comparison_tsv), delimiter="\t"))
    if not rows:
        raise SystemExit("No rows found")

    by_block = defaultdict(list)

    for i, r in enumerate(rows):
        pred = ival(r.get("pred_state", r.get("local_phase_state", 0)), 0)
        truth = ival(r.get("truth_state", 0), 0)
        block = str(r.get("block_id", "unassigned"))
        if pred != 0 and truth != 0 and block not in {"", "unassigned", "0"}:
            by_block[block].append((i, pred, truth))

    orientation = {}
    counts = {}

    for block, vals in by_block.items():
        same = sum(1 for _, p, t in vals if p == t)
        opposite = sum(1 for _, p, t in vals if -p == t)

        if same > opposite:
            orient = 1
        elif opposite > same:
            orient = -1
        else:
            orient = 0

        orientation[block] = orient
        counts[block] = (same, opposite, abs(same - opposite))

    out_rows = []
    correct = wrong = not_comparable = unknown_orientation = 0

    for r in rows:
        rr = dict(r)

        pred = ival(rr.get("pred_state", rr.get("local_phase_state", 0)), 0)
        truth = ival(rr.get("truth_state", 0), 0)
        block = str(rr.get("block_id", "unassigned"))

        orient = orientation.get(block, 0)
        same, opposite, margin = counts.get(block, (0, 0, 0))

        oriented_pred = pred * orient if pred != 0 and orient != 0 else 0

        comparable = pred != 0 and truth != 0
        if not comparable:
            corr = ""
            not_comparable += 1
        elif orient == 0:
            corr = ""
            unknown_orientation += 1
        else:
            corr = 1 if oriented_pred == truth else 0
            if corr:
                correct += 1
            else:
                wrong += 1

        rr["block_orientation"] = orient
        rr["oriented_pred_state"] = oriented_pred
        rr["correct_after_block_orientation"] = corr
        rr["orientation_same_count"] = same
        rr["orientation_opposite_count"] = opposite
        rr["orientation_margin"] = margin

        out_rows.append(rr)

    fields = list(rows[0].keys())
    for f in [
        "block_orientation",
        "oriented_pred_state",
        "correct_after_block_orientation",
        "orientation_same_count",
        "orientation_opposite_count",
        "orientation_margin",
    ]:
        if f not in fields:
            fields.append(f)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    denom = correct + wrong
    acc = correct / denom if denom else ""

    with open(args.out_summary_tsv, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"correct\t{correct}\n")
        fh.write(f"wrong\t{wrong}\n")
        fh.write(f"not_comparable\t{not_comparable}\n")
        fh.write(f"unknown_orientation\t{unknown_orientation}\n")
        fh.write(f"oriented_site_accuracy\t{acc}\n")

    print({
        "rows": len(rows),
        "blocks_oriented": len([b for b, o in orientation.items() if o != 0]),
        "correct": correct,
        "wrong": wrong,
        "unknown_orientation": unknown_orientation,
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
