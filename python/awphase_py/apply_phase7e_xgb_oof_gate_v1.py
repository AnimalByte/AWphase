#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

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

def is_panel_filled(r):
    return str(r.get("phase7a_panel_filled", "0")).strip().lower() in {"1", "true", "t", "yes"}

def zero_call(r):
    if "local_phase_state" in r:
        r["local_phase_state"] = 0
    if "phase_state" in r:
        r["phase_state"] = 0
    r["block_id"] = "unassigned"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--xgb-predictions-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    pred = {}
    for r in read_tsv(args.xgb_predictions_tsv):
        if r.get("source_window") != args.source_window:
            continue
        try:
            pos = int(float(r["pos"]))
        except Exception:
            continue
        pred[pos] = fval(r.get("xgb_oof_probability"), 0.0)

    rows = read_tsv(args.local_calls_tsv)

    total_panel_filled = 0
    kept = 0
    removed = 0
    no_xgb_prediction = 0

    out = []
    for r in rows:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)

        if pos is not None and is_panel_filled(rr):
            total_panel_filled += 1
            p = pred.get(pos)

            if p is None:
                no_xgb_prediction += 1
                zero_call(rr)
                rr["phase7e_xgb_kept"] = 0
                rr["phase7e_xgb_probability"] = ""
                rr["phase7e_xgb_reason"] = "no_oof_prediction"
                removed += 1
            elif p >= args.threshold:
                kept += 1
                rr["phase7e_xgb_kept"] = 1
                rr["phase7e_xgb_probability"] = f"{p:.6f}"
                rr["phase7e_xgb_reason"] = "kept"
            else:
                removed += 1
                zero_call(rr)
                rr["phase7e_xgb_kept"] = 0
                rr["phase7e_xgb_probability"] = f"{p:.6f}"
                rr["phase7e_xgb_reason"] = "below_threshold"
        else:
            rr["phase7e_xgb_kept"] = ""
            rr["phase7e_xgb_probability"] = ""
            rr["phase7e_xgb_reason"] = ""

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
        "source_window": args.source_window,
        "threshold": args.threshold,
        "total_panel_filled": total_panel_filled,
        "kept": kept,
        "removed": removed,
        "no_xgb_prediction": no_xgb_prediction,
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
