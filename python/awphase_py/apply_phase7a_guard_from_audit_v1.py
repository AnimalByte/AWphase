#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

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

def is_panel_filled(r):
    return str(r.get("phase7a_panel_filled", "0")).strip().lower() in {"1", "true", "t", "yes"}

def set_zero(r):
    if "local_phase_state" in r:
        r["local_phase_state"] = 0
    if "phase_state" in r:
        r["phase_state"] = 0
    r["block_id"] = "unassigned"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--audit-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)

    ap.add_argument("--min-orientation-sites", type=int, default=0)
    ap.add_argument("--min-orientation-margin", type=int, default=0)
    ap.add_argument("--min-panel-confidence", type=float, default=0.0)
    ap.add_argument("--min-panel-margin", type=float, default=0.0)
    ap.add_argument("--min-best-vs-second-margin", type=float, default=0.0)
    ap.add_argument("--min-anchors", type=int, default=0)
    args = ap.parse_args()

    audit = {}
    for r in read_tsv(args.audit_tsv):
        try:
            pos = int(float(r["pos"]))
        except Exception:
            continue
        audit[pos] = r

    rows = read_tsv(args.local_calls_tsv)

    kept = 0
    removed = 0
    total_filled = 0

    out = []
    for r in rows:
        rr = dict(r)
        try:
            pos = int(float(rr["pos"]))
        except Exception:
            out.append(rr)
            continue

        if is_panel_filled(rr):
            total_filled += 1
            ar = audit.get(pos, {})

            ok = True
            if ival(ar.get("orientation_sites"), 0) < args.min_orientation_sites:
                ok = False
            if ival(ar.get("orientation_margin"), 0) < args.min_orientation_margin:
                ok = False
            if fval(ar.get("panel_confidence"), 0.0) < args.min_panel_confidence:
                ok = False
            if fval(ar.get("panel_margin"), 0.0) < args.min_panel_margin:
                ok = False
            if fval(ar.get("best_vs_second_margin"), 0.0) < args.min_best_vs_second_margin:
                ok = False
            if ival(ar.get("anchors"), 0) < args.min_anchors:
                ok = False

            if ok:
                kept += 1
                rr["phase7a_guard_kept"] = 1
            else:
                removed += 1
                set_zero(rr)
                rr["phase7a_guard_kept"] = 0
                rr["phase7a_guard_removed"] = 1
        else:
            rr["phase7a_guard_kept"] = ""
            rr["phase7a_guard_removed"] = ""

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
        "total_panel_filled": total_filled,
        "kept": kept,
        "removed": removed,
        "min_orientation_sites": args.min_orientation_sites,
        "min_orientation_margin": args.min_orientation_margin,
        "min_panel_confidence": args.min_panel_confidence,
        "min_panel_margin": args.min_panel_margin,
        "min_best_vs_second_margin": args.min_best_vs_second_margin,
        "min_anchors": args.min_anchors,
    })

if __name__ == "__main__":
    main()
