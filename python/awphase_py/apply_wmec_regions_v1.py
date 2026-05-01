#!/usr/bin/env python3
import argparse, csv
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    seen = set(fields); final = list(fields)
    for r in rows:
        for k in r:
            if k not in seen:
                final.append(k); seen.add(k)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=final, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def as_int(x, d=0):
    try: return int(float(x))
    except Exception: return d


def as_float(x, d=0.0):
    try: return float(x)
    except Exception: return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--wmec-solutions-tsv", required=True)
    ap.add_argument("--apply-margin", type=float, default=0.55)
    ap.add_argument("--min-real-fragments", type=int, default=1)
    ap.add_argument("--min-total-fragments", type=int, default=1)
    ap.add_argument("--require-anchor-or-donor", action='store_true')
    ap.add_argument("--min-abs-donor", type=float, default=0.03)
    ap.add_argument("--exact-only", action='store_true')
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    sol_by_pos = {}
    for r in read_tsv(args.wmec_solutions_tsv):
        pos = as_int(r.get("pos"))
        if not pos: continue
        real = as_int(r.get("real_fragments", 0))
        scaf = as_int(r.get("scaffold_fragments", 0))
        sol_by_pos[pos] = {
            "state": as_int(r.get("resolved_state", 0)),
            "margin": as_float(r.get("region_margin", 0.0)),
            "anchor_ok": as_int(r.get("anchor_ok", 0)),
            "avg_donor_bias": as_float(r.get("avg_donor_bias", 0.0)),
            "region_id": r.get("region_id", ""),
            "used_exact": as_int(r.get("used_exact", 0)),
            "real_fragments": real,
            "scaffold_fragments": scaf,
            "total_fragments": real + scaf,
        }

    rows = read_tsv(args.local_calls_tsv)
    out = []; modified = 0; considered = 0; gated = {"low_margin":0,"low_fragments":0,"anchor_donor":0,"exact":0}
    for r in rows:
        pos = as_int(r.get("pos"))
        row = dict(r)
        for k, v in {
            "wmec_applied": 0,
            "wmec_region_id": "",
            "wmec_margin": "",
            "wmec_anchor_ok": "",
            "wmec_avg_donor_bias": "",
            "wmec_used_exact": "",
            "wmec_real_fragments": "",
            "wmec_scaffold_fragments": "",
            "wmec_gate_reason": "",
        }.items():
            row[k] = v
        s = sol_by_pos.get(pos)
        if s:
            considered += 1
            reasons = []
            if s["margin"] < args.apply_margin:
                reasons.append("low_margin"); gated["low_margin"] += 1
            if s["real_fragments"] < args.min_real_fragments or s["total_fragments"] < args.min_total_fragments:
                reasons.append("low_fragments"); gated["low_fragments"] += 1
            if args.exact_only and s["used_exact"] != 1:
                reasons.append("not_exact"); gated["exact"] += 1
            donor_ok = abs(s["avg_donor_bias"]) >= args.min_abs_donor
            if args.require_anchor_or_donor and not (s["anchor_ok"] == 1 or donor_ok):
                reasons.append("anchor_donor"); gated["anchor_donor"] += 1
            row["wmec_region_id"] = s["region_id"]
            row["wmec_margin"] = s["margin"]
            row["wmec_anchor_ok"] = s["anchor_ok"]
            row["wmec_avg_donor_bias"] = s["avg_donor_bias"]
            row["wmec_used_exact"] = s["used_exact"]
            row["wmec_real_fragments"] = s["real_fragments"]
            row["wmec_scaffold_fragments"] = s["scaffold_fragments"]
            if not reasons:
                target_col = "local_phase_state" if "local_phase_state" in row else "phase_state"
                old = as_int(row.get(target_col, 0))
                if old != s["state"]:
                    row[target_col] = s["state"]
                    row["wmec_applied"] = 1
                    modified += 1
                row["wmec_gate_reason"] = "applied_or_same"
            else:
                row["wmec_gate_reason"] = "|".join(reasons)
        out.append(row)

    write_tsv(args.out_tsv, out, list(rows[0].keys()) if rows else [])
    print({"rows": len(out), "considered_positions": considered, "modified_rows": modified, "apply_margin": args.apply_margin, "gated": gated})

if __name__ == "__main__":
    main()
