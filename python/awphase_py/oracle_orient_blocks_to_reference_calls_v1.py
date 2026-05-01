#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                fields.append(k)
                seen.add(k)

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def write_report(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "block_id",
        "n_target_nonzero",
        "n_ref_overlap",
        "same_count",
        "opposite_count",
        "chosen_orientation",
        "margin",
        "flipped",
        "reason",
    ]

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-tsv", required=True, help="AWPhase local_calls.tsv to re-orient diagnostically")
    ap.add_argument("--reference-tsv", required=True, help="Reference local_calls.tsv, e.g. WhatsHap")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-report-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    ap.add_argument("--min-overlap-sites", type=int, default=2)
    ap.add_argument("--min-margin", type=int, default=1)
    args = ap.parse_args()

    target = read_tsv(args.target_tsv)
    ref_rows = read_tsv(args.reference_tsv)

    ref_state = {}
    for r in ref_rows:
        try:
            p = int(float(r["pos"]))
        except Exception:
            continue
        s = state(r)
        if s != 0:
            ref_state[p] = s

    by_block = defaultdict(list)
    for i, r in enumerate(target):
        try:
            p = int(float(r["pos"]))
        except Exception:
            continue
        s = state(r)
        b = block_id(r)
        if s != 0 and b:
            by_block[b].append((i, p, s))

    orientation = {}
    report = []

    blocks_considered = 0
    blocks_flipped = 0
    blocks_kept = 0
    blocks_unoriented = 0
    overlap_sites = 0
    same_total = 0
    opposite_total = 0

    for b, items in sorted(by_block.items()):
        blocks_considered += 1
        same = 0
        opp = 0
        n_overlap = 0

        for _, p, s in items:
            rs = ref_state.get(p)
            if rs is None:
                continue
            n_overlap += 1
            if s == rs:
                same += 1
            elif s == -rs:
                opp += 1

        overlap_sites += n_overlap
        same_total += same
        opposite_total += opp

        margin = abs(same - opp)
        chosen = 1
        reason = "keep_default"

        if n_overlap < args.min_overlap_sites:
            chosen = 1
            reason = "too_few_overlap"
            blocks_unoriented += 1
        elif margin < args.min_margin:
            chosen = 1
            reason = "low_margin"
            blocks_unoriented += 1
        elif opp > same:
            chosen = -1
            reason = "oracle_flip"
            blocks_flipped += 1
        else:
            chosen = 1
            reason = "oracle_keep"
            blocks_kept += 1

        orientation[b] = chosen

        report.append({
            "block_id": b,
            "n_target_nonzero": len(items),
            "n_ref_overlap": n_overlap,
            "same_count": same,
            "opposite_count": opp,
            "chosen_orientation": chosen,
            "margin": margin,
            "flipped": int(chosen == -1),
            "reason": reason,
        })

    out = [dict(r) for r in target]
    changed_sites = 0

    for b, items in by_block.items():
        o = orientation.get(b, 1)
        if o == 1:
            continue
        for i, p, s in items:
            set_state(out[i], -s)
            out[i]["oracle_reference_orientation"] = -1
            changed_sites += 1

    for r in out:
        r.setdefault("oracle_reference_orientation", "1")

    write_tsv(args.out_tsv, out)
    write_report(args.out_report_tsv, report)

    summary = {
        "target_tsv": args.target_tsv,
        "reference_tsv": args.reference_tsv,
        "blocks_considered": blocks_considered,
        "blocks_flipped": blocks_flipped,
        "blocks_kept": blocks_kept,
        "blocks_unoriented": blocks_unoriented,
        "overlap_sites_with_reference": overlap_sites,
        "same_total_before": same_total,
        "opposite_total_before": opposite_total,
        "same_fraction_before": same_total / overlap_sites if overlap_sites else None,
        "opposite_fraction_before": opposite_total / overlap_sites if overlap_sites else None,
        "changed_sites": changed_sites,
        "min_overlap_sites": args.min_overlap_sites,
        "min_margin": args.min_margin,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
