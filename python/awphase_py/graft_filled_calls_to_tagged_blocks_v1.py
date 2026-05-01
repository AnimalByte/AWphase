#!/usr/bin/env python3
import argparse
import bisect
import csv
import json
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

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def fval(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def state(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def set_state(row, s):
    row["local_phase_state"] = int(s)
    row["phase_state"] = int(s)

def usable_block_id(row):
    bid = str(row.get("block_id", "")).strip()
    if bid == "" or bid.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return bid

def nearest_anchors(anchor_positions, pos, max_dist):
    i = bisect.bisect_left(anchor_positions, pos)
    left = None
    right = None

    if i > 0:
        lp = anchor_positions[i - 1]
        if pos - lp <= max_dist:
            left = lp

    if i < len(anchor_positions):
        rp = anchor_positions[i]
        if rp - pos <= max_dist:
            right = rp

    return left, right

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tagged-local-calls-tsv", required=True)
    ap.add_argument("--filled-local-calls-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--max-anchor-dist-bp", type=int, default=100000)
    ap.add_argument(
        "--mode",
        choices=["same_block_only", "nearest_block", "left_prefer"],
        default="same_block_only",
    )
    ap.add_argument("--require-state-agreement-with-anchor", action="store_true")
    args = ap.parse_args()

    tagged = read_tsv(args.tagged_local_calls_tsv)
    filled = read_tsv(args.filled_local_calls_tsv)

    tagged_by_pos = {ival(r["pos"]): r for r in tagged}
    filled_by_pos = {ival(r["pos"]): r for r in filled}

    # Backbone anchors are tagged nonzero calls with real block IDs.
    anchor_state = {}
    anchor_block = {}
    for pos, r in tagged_by_pos.items():
        st = state(r)
        bid = usable_block_id(r)
        if st != 0 and bid:
            anchor_state[pos] = st
            anchor_block[pos] = bid

    anchor_positions = sorted(anchor_state)

    out = [dict(r) for r in filled]
    out_by_pos = {ival(r["pos"]): r for r in out}

    considered = 0
    grafted = 0
    changed_block = 0
    skipped_no_tagged_row = 0
    skipped_not_fill = 0
    skipped_no_anchor = 0
    skipped_cross_block = 0
    skipped_state_guard = 0

    for pos, row in out_by_pos.items():
        tag_row = tagged_by_pos.get(pos)
        if tag_row is None:
            skipped_no_tagged_row += 1
            continue

        tag_st = state(tag_row)
        fill_st = state(row)

        # Only graft newly filled sites: tagged was abstain, filled has state.
        if tag_st != 0 or fill_st == 0:
            skipped_not_fill += 1
            continue

        considered += 1

        left, right = nearest_anchors(anchor_positions, pos, args.max_anchor_dist_bp)
        if left is None and right is None:
            skipped_no_anchor += 1
            continue

        chosen_block = ""
        chosen_anchor = None

        if args.mode == "same_block_only":
            if left is not None and right is not None and anchor_block[left] == anchor_block[right]:
                chosen_block = anchor_block[left]
                chosen_anchor = left if abs(pos - left) <= abs(right - pos) else right
            elif left is not None and right is None:
                chosen_block = anchor_block[left]
                chosen_anchor = left
            elif right is not None and left is None:
                chosen_block = anchor_block[right]
                chosen_anchor = right
            else:
                skipped_cross_block += 1
                continue

        elif args.mode == "nearest_block":
            if left is not None and right is not None:
                chosen_anchor = left if abs(pos - left) <= abs(right - pos) else right
            else:
                chosen_anchor = left if left is not None else right
            chosen_block = anchor_block[chosen_anchor]

        elif args.mode == "left_prefer":
            chosen_anchor = left if left is not None else right
            chosen_block = anchor_block[chosen_anchor]

        if args.require_state_agreement_with_anchor:
            # This is conservative. It only grafts if the filled state agrees with the chosen tagged anchor state.
            # If the fill source is already read-bridge oriented, this should usually pass for high-confidence intervals.
            if chosen_anchor is not None and fill_st != anchor_state[chosen_anchor]:
                skipped_state_guard += 1
                continue

        old_block = str(row.get("block_id", "")).strip()
        row["phase5b_old_block_id"] = old_block
        row["phase5b_grafted_block_id"] = chosen_block
        row["phase5b_anchor_pos"] = "" if chosen_anchor is None else chosen_anchor
        row["phase5b_anchor_state"] = "" if chosen_anchor is None else anchor_state[chosen_anchor]
        row["phase5b_mode"] = args.mode

        row["block_id"] = chosen_block
        set_state(row, fill_st)

        grafted += 1
        if old_block != chosen_block:
            changed_block += 1

    summary = {
        "rows": len(out),
        "backbone_anchors": len(anchor_positions),
        "considered_new_fill_sites": considered,
        "grafted_sites": grafted,
        "changed_block": changed_block,
        "skipped_no_tagged_row": skipped_no_tagged_row,
        "skipped_not_fill": skipped_not_fill,
        "skipped_no_anchor": skipped_no_anchor,
        "skipped_cross_block": skipped_cross_block,
        "skipped_state_guard": skipped_state_guard,
        "max_anchor_dist_bp": args.max_anchor_dist_bp,
        "mode": args.mode,
        "require_state_agreement_with_anchor": bool(args.require_state_agreement_with_anchor),
    }

    write_tsv(args.out_tsv, out)
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
