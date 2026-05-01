#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def load_site_comparison(path):
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = None
            for key in ("pos", "site_pos", "position", "variant_pos"):
                if key in row and str(row[key]).strip():
                    pos = int(str(row[key]).strip())
                    break
            if pos is None:
                continue
            row["_pos"] = pos
            rows.append(row)
    return rows

def find_sites_payload(obj):
    if isinstance(obj, dict) and "sites" in obj and isinstance(obj["sites"], list):
        return obj["sites"]
    raise SystemExit("Could not find 'sites' list in local debug JSON")

def row_wrong(row):
    for key in ("is_correct", "correct", "pred_matches_truth", "match"):
        if key in row:
            val = str(row[key]).strip().lower()
            if val in {"false", "0", "no", "wrong", "incorrect", "mismatch"}:
                return True
            if val in {"true", "1", "yes", "correct"}:
                return False
    pred = None
    truth = None
    for key in ("pred_phase_state", "predicted_phase_state", "phase_state_pred", "pred_state"):
        if key in row and str(row[key]).strip():
            pred = str(row[key]).strip()
            break
    for key in ("truth_phase_state", "phase_state_truth", "truth_state"):
        if key in row and str(row[key]).strip():
            truth = str(row[key]).strip()
            break
    return pred is not None and truth is not None and pred != truth

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-debug", required=True)
    ap.add_argument("--site-comparison", required=True)
    ap.add_argument("--window-index", type=int, action="append", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    debug = load_json(args.local_debug)
    sites = find_sites_payload(debug)
    cmp_rows = load_site_comparison(args.site_comparison)
    cmp_by_pos = {r["_pos"]: r for r in cmp_rows}

    selected = []
    wanted = set(args.window_index)
    for s in sites:
        wi = int(s.get("window_index", -1))
        if wi not in wanted:
            continue
        pos = int(s["pos"])
        cmp_row = cmp_by_pos.get(pos, {})
        selected.append({
            "window_index": wi,
            "pos": pos,
            "phase_state": s.get("phase_state"),
            "confidence": s.get("confidence"),
            "donor_margin": s.get("donor_margin"),
            "path_sign": s.get("path_sign"),
            "path_margin_agrees": s.get("path_margin_agrees"),
            "decision_source": s.get("decision_source"),
            "truth_row_present": bool(cmp_row),
            "wrong_vs_truth": row_wrong(cmp_row) if cmp_row else "",
            "pred_phase_state": cmp_row.get("pred_phase_state", cmp_row.get("predicted_phase_state", "")),
            "truth_phase_state": cmp_row.get("truth_phase_state", cmp_row.get("truth_state", "")),
        })

    selected.sort(key=lambda r: (r["window_index"], r["pos"]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "window_index",
                "pos",
                "phase_state",
                "confidence",
                "donor_margin",
                "path_sign",
                "path_margin_agrees",
                "decision_source",
                "truth_row_present",
                "wrong_vs_truth",
                "pred_phase_state",
                "truth_phase_state",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(selected)

    print(f"Wrote {len(selected)} rows to {out}")

if __name__ == "__main__":
    main()
