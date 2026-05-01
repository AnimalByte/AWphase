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
    ap.add_argument("--out-tsv", required=True)

    ap.add_argument("--min-region-margin", type=float, default=0.05)
    ap.add_argument("--min-orientation-margin", type=float, default=0.75)
    ap.add_argument("--min-donor-abs", type=float, default=0.05)
    ap.add_argument("--require-anchor-or-donor", action="store_true")
    ap.add_argument("--allow-singletons", action="store_true")
    args = ap.parse_args()

    local_rows = read_tsv(args.local_calls_tsv)
    by_pos = {ival(r["pos"]): r for r in local_rows}

    sols = read_tsv(args.oriented_resolutions_tsv)

    region_sizes = {}
    for r in sols:
        rid = r.get("region_id", "")
        region_sizes[rid] = region_sizes.get(rid, 0) + 1

    considered = 0
    applied_rows = 0
    changed_rows = 0
    skipped_margin = 0
    skipped_orientation = 0
    skipped_singleton = 0
    skipped_guard = 0
    missing_pos = 0

    for s in sols:
        rid = s.get("region_id", "")
        pos = ival(s.get("pos"))
        state = ival(s.get("resolved_state"), 0)

        if pos not in by_pos:
            missing_pos += 1
            continue
        if state == 0:
            continue

        considered += 1

        region_margin = fval(s.get("region_margin"), 0.0)
        orientation_margin = fval(s.get("orientation_margin"), 0.0)
        donor = fval(s.get("avg_donor_bias"), 0.0)
        anchor_ok = ival(s.get("anchor_ok"), 0)

        if region_margin < args.min_region_margin:
            skipped_margin += 1
            continue

        if orientation_margin < args.min_orientation_margin:
            skipped_orientation += 1
            continue

        if not args.allow_singletons and region_sizes.get(rid, 0) <= 1:
            skipped_singleton += 1
            continue

        if args.require_anchor_or_donor and not (anchor_ok == 1 or abs(donor) >= args.min_donor_abs):
            skipped_guard += 1
            continue

        row = by_pos[pos]
        old = ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

        row["phase4d_old_local_phase_state"] = old
        row["phase4d_region_id"] = rid
        row["phase4d_resolved_state"] = state
        row["phase4d_resolved_state_original"] = s.get("resolved_state_original", "")
        row["phase4d_region_orientation"] = s.get("region_orientation", "")
        row["phase4d_region_margin"] = f"{region_margin:.6f}"
        row["phase4d_orientation_margin"] = f"{orientation_margin:.6f}"
        row["phase4d_anchor_ok"] = anchor_ok
        row["phase4d_avg_donor_bias"] = f"{donor:.6f}"
        row["phase4d_applied"] = 1

        # Critical: evaluator reads local_phase_state.
        row["local_phase_state"] = state
        row["phase_state"] = state

        applied_rows += 1
        if old != state:
            changed_rows += 1

    for r in local_rows:
        r.setdefault("phase4d_applied", 0)

    write_tsv(args.out_tsv, local_rows)

    print({
        "rows": len(local_rows),
        "solutions_considered": considered,
        "applied_rows": applied_rows,
        "changed_rows": changed_rows,
        "missing_pos": missing_pos,
        "skipped_margin": skipped_margin,
        "skipped_orientation": skipped_orientation,
        "skipped_singleton": skipped_singleton,
        "skipped_guard": skipped_guard,
        "min_region_margin": args.min_region_margin,
        "min_orientation_margin": args.min_orientation_margin,
        "min_donor_abs": args.min_donor_abs,
        "require_anchor_or_donor": bool(args.require_anchor_or_donor),
        "allow_singletons": bool(args.allow_singletons),
    })

if __name__ == "__main__":
    main()
