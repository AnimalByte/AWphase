#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict, Counter

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-comparison-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-tsv", required=True)
    ap.add_argument("--min-orientation-sites", type=int, default=2)
    args = ap.parse_args()

    rows = read_tsv(args.site_comparison_tsv)

    counts = defaultdict(lambda: {"same": 0, "opposite": 0})

    for r in rows:
        block = str(r.get("block_id", "")).strip()
        pred = ival(r.get("pred_state"), 0)
        truth = ival(r.get("truth_state"), 0)
        matched = str(r.get("matched_exact_truth", "1")).strip().lower()

        if not block or block.lower() in {"unassigned", "none", "na", "."}:
            continue
        if matched in {"0", "false"}:
            continue
        if pred == 0 or truth == 0:
            continue

        if pred == truth:
            counts[block]["same"] += 1
        elif pred == -truth:
            counts[block]["opposite"] += 1

    orientation = {}
    for block, c in counts.items():
        same = c["same"]
        opp = c["opposite"]
        n = same + opp

        if n < args.min_orientation_sites:
            orientation[block] = ""
        else:
            orientation[block] = 1 if same >= opp else -1

    out = []
    summary = Counter()

    for r in rows:
        rr = dict(r)

        block = str(r.get("block_id", "")).strip()
        pred = ival(r.get("pred_state"), 0)
        truth = ival(r.get("truth_state"), 0)
        matched = str(r.get("matched_exact_truth", "1")).strip().lower()

        o = orientation.get(block, "")
        oriented_pred = ""
        correct = ""

        if pred == 0 or truth == 0 or matched in {"0", "false"}:
            correct = ""
            summary["not_comparable"] += 1
        elif o == "":
            correct = ""
            summary["unknown_orientation"] += 1
        else:
            oriented_pred = pred * int(o)
            correct = int(oriented_pred == truth)
            if correct:
                summary["correct"] += 1
            else:
                summary["wrong"] += 1

        rr["block_orientation"] = o
        rr["oriented_pred_state"] = oriented_pred
        rr["correct_after_block_orientation"] = correct
        rr["orientation_same_count"] = counts.get(block, {}).get("same", 0)
        rr["orientation_opposite_count"] = counts.get(block, {}).get("opposite", 0)
        rr["orientation_margin"] = abs(
            counts.get(block, {}).get("same", 0)
            - counts.get(block, {}).get("opposite", 0)
        )

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

    with open(args.out_summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"], delimiter="\t")
        w.writeheader()
        for k in ["correct", "wrong", "not_comparable", "unknown_orientation"]:
            w.writerow({"metric": k, "value": summary.get(k, 0)})

        denom = summary.get("correct", 0) + summary.get("wrong", 0)
        w.writerow({
            "metric": "oriented_site_accuracy",
            "value": summary.get("correct", 0) / denom if denom else "",
        })

    print(dict(summary))

if __name__ == "__main__":
    main()
