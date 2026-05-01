#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

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

def fval(x, default=0.0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return float(x)
    except Exception:
        return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--xgb-predictions-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    base = read_tsv(args.local_calls_phase6c_tsv)
    preds = read_tsv(args.xgb_predictions_tsv)

    best_by_pos = {}

    for r in preds:
        if str(r.get("source_window")) != args.source_window:
            continue

        p = fval(r.get("xgb_oof_probability"), -1.0)
        if p < args.threshold:
            continue

        pos = ival(r.get("pos"), None)
        pred_state = ival(r.get("pred_state"), 0)
        block = str(r.get("selected_block_id", "unassigned"))

        if pos is None or pred_state == 0 or block in {"", "0", "unassigned"}:
            continue

        old = best_by_pos.get(pos)
        if old is None or p > old["prob"]:
            best_by_pos[pos] = {
                "prob": p,
                "pred_state": pred_state,
                "selected_block_id": block,
            }

    out_rows = []
    projected = 0

    for r in base:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        state = ival(rr.get("local_phase_state", rr.get("phase_state", 0)), 0)

        rr["phase8d_projected"] = 0
        rr["phase8d_xgb_probability"] = ""

        if pos in best_by_pos and state == 0:
            hit = best_by_pos[pos]
            if "local_phase_state" in rr:
                rr["local_phase_state"] = hit["pred_state"]
            if "phase_state" in rr:
                rr["phase_state"] = hit["pred_state"]
            rr["block_id"] = hit["selected_block_id"]
            rr["phase8d_projected"] = 1
            rr["phase8d_xgb_probability"] = hit["prob"]
            projected += 1

        out_rows.append(rr)

    fields = list(base[0].keys()) if base else []
    for f in ["phase8d_projected", "phase8d_xgb_probability"]:
        if f not in fields:
            fields.append(f)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    print({
        "source_window": args.source_window,
        "threshold": args.threshold,
        "accepted_projected_sites": projected,
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
