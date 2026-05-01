#!/usr/bin/env python3
import argparse
import csv
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--oriented-resolutions-tsv", required=True)
    ap.add_argument("--region-orientations-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)

    ap.add_argument("--min-region-margin", type=float, default=0.05)
    ap.add_argument("--min-flip-orientation-margin", type=float, default=1.0)
    ap.add_argument("--min-bridge-read-pairs-for-flip", type=int, default=2)
    ap.add_argument("--min-real-fragments", type=int, default=1)
    ap.add_argument("--allow-singletons", action="store_true")
    args = ap.parse_args()

    local_rows = read_tsv(args.local_calls_tsv)
    by_pos = {ival(r["pos"]): r for r in local_rows}

    orient_rows = read_tsv(args.region_orientations_tsv)
    orient = {}
    for r in orient_rows:
        rid = r.get("region_id", "")
        if not rid:
            continue
        orient[rid] = {
            "orientation": ival(r.get("orientation"), 1),
            "orientation_margin": fval(r.get("orientation_margin"), 0.0),
            "bridge_read_pairs": ival(r.get("bridge_read_pairs"), 0),
            "n_sites": ival(r.get("n_sites"), 0),
        }

    sols = read_tsv(args.oriented_resolutions_tsv)

    region_sizes = {}
    for r in sols:
        rid = r.get("region_id", "")
        region_sizes[rid] = region_sizes.get(rid, 0) + 1

    considered = 0
    applied = 0
    changed = 0
    proposed_flips = 0
    accepted_flips = 0
    rejected_flips_no_bridge = 0
    rejected_flips_low_margin = 0
    skipped_margin = 0
    skipped_singleton = 0
    skipped_real_fragments = 0

    for s in sols:
        rid = s.get("region_id", "")
        pos = ival(s.get("pos"))
        if pos not in by_pos:
            continue

        region_margin = fval(s.get("region_margin"), 0.0)
        if region_margin < args.min_region_margin:
            skipped_margin += 1
            continue

        if not args.allow_singletons and region_sizes.get(rid, 0) <= 1:
            skipped_singleton += 1
            continue

        real_fragments = ival(s.get("real_fragments"), 0)
        if real_fragments < args.min_real_fragments:
            skipped_real_fragments += 1
            continue

        oriented_state = ival(s.get("resolved_state"), 0)
        original_state = ival(s.get("resolved_state_original", oriented_state), oriented_state)
        if original_state == 0 and oriented_state != 0:
            original_state = oriented_state

        if oriented_state == 0:
            continue

        o = orient.get(rid, {})
        orientation = o.get("orientation", 1)
        orientation_margin = o.get("orientation_margin", fval(s.get("orientation_margin"), 0.0))
        bridge_pairs = o.get("bridge_read_pairs", 0)

        final_state = original_state

        if orientation == -1:
            proposed_flips += 1
            if bridge_pairs < args.min_bridge_read_pairs_for_flip:
                rejected_flips_no_bridge += 1
            elif orientation_margin < args.min_flip_orientation_margin:
                rejected_flips_low_margin += 1
            else:
                final_state = -original_state
                accepted_flips += 1
        else:
            final_state = original_state

        row = by_pos[pos]
        old = ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

        row["phase4e_old_local_phase_state"] = old
        row["phase4e_region_id"] = rid
        row["phase4e_original_state"] = original_state
        row["phase4e_oriented_state"] = oriented_state
        row["phase4e_final_state"] = final_state
        row["phase4e_region_orientation"] = orientation
        row["phase4e_orientation_margin"] = f"{orientation_margin:.6f}"
        row["phase4e_bridge_read_pairs"] = bridge_pairs
        row["phase4e_region_margin"] = f"{region_margin:.6f}"
        row["phase4e_real_fragments"] = real_fragments
        row["phase4e_applied"] = 1

        row["local_phase_state"] = final_state
        row["phase_state"] = final_state

        considered += 1
        applied += 1
        if old != final_state:
            changed += 1

    for r in local_rows:
        r.setdefault("phase4e_applied", 0)

    write_tsv(args.out_tsv, local_rows)

    print({
        "rows": len(local_rows),
        "solutions_considered": considered,
        "applied_rows": applied,
        "changed_rows": changed,
        "proposed_flip_rows": proposed_flips,
        "accepted_flip_rows": accepted_flips,
        "rejected_flips_no_bridge": rejected_flips_no_bridge,
        "rejected_flips_low_margin": rejected_flips_low_margin,
        "skipped_margin": skipped_margin,
        "skipped_singleton": skipped_singleton,
        "skipped_real_fragments": skipped_real_fragments,
        "min_region_margin": args.min_region_margin,
        "min_flip_orientation_margin": args.min_flip_orientation_margin,
        "min_bridge_read_pairs_for_flip": args.min_bridge_read_pairs_for_flip,
        "min_real_fragments": args.min_real_fragments,
        "allow_singletons": bool(args.allow_singletons),
    })

if __name__ == "__main__":
    main()
