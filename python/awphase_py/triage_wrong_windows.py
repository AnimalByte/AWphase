#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def walk(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from walk(item)

def find_debug_sites(obj):
    if isinstance(obj, dict) and isinstance(obj.get("sites"), list):
        return obj["sites"]
    for d in walk(obj):
        if isinstance(d, dict) and isinstance(d.get("sites"), list):
            return d["sites"]
    return []

def parse_pos(row):
    for key in ("pos", "site_pos", "position", "variant_pos"):
        if key in row and str(row[key]).strip():
            try:
                return int(str(row[key]).strip())
            except ValueError:
                pass
    return None

def get_pred_state(row):
    for key in ("pred_phase_state", "predicted_phase_state", "phase_state_pred", "pred_state", "local_phase_state"):
        if key in row and str(row[key]).strip():
            try:
                return int(str(row[key]).strip())
            except ValueError:
                return str(row[key]).strip()
    return None

def get_truth_state(row):
    for key in ("truth_phase_state", "phase_state_truth", "truth_state"):
        if key in row and str(row[key]).strip():
            try:
                return int(str(row[key]).strip())
            except ValueError:
                return str(row[key]).strip()
    return None

def classify_row(row):
    pred = get_pred_state(row)
    truth = get_truth_state(row)

    if pred is None or truth is None:
        return "unknown"

    if pred == 0:
        return "abstain"

    if pred == truth:
        return "correct"

    return "wrong_nonzero"

def mean(values):
    return sum(values) / len(values) if values else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-debug", required=True)
    ap.add_argument("--site-comparison", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    debug_obj = load_json(args.local_debug)
    sites = find_debug_sites(debug_obj)

    pos_class = {}
    with open(args.site_comparison) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = parse_pos(row)
            if pos is None:
                continue
            pos_class[pos] = classify_row(row)

    by_window = defaultdict(list)
    for s in sites:
        try:
            w = int(s.get("window_index", 0))
        except Exception:
            w = 0
        by_window[w].append(s)

    rows = []
    for w in sorted(by_window):
        recs = by_window[w]
        positions = [int(r["pos"]) for r in recs if "pos" in r and str(r["pos"]).strip()]

        wrong_nonzero = [p for p in positions if pos_class.get(p) == "wrong_nonzero"]
        abstain = [p for p in positions if pos_class.get(p) == "abstain"]
        correct = [p for p in positions if pos_class.get(p) == "correct"]

        confidences = [float(r["confidence"]) for r in recs if str(r.get("confidence", "")).strip() not in ("", "None")]
        donor_margins = [abs(float(r["donor_margin"])) for r in recs if str(r.get("donor_margin", "")).strip() not in ("", "None")]
        decision_sources = sorted({str(r.get("decision_source", "")) for r in recs if str(r.get("decision_source", "")).strip()})

        rows.append({
            "window_index": w,
            "n_sites": len(recs),
            "n_wrong_nonzero": len(wrong_nonzero),
            "n_abstain": len(abstain),
            "n_correct": len(correct),
            "wrong_nonzero_fraction": f"{(len(wrong_nonzero) / len(recs)) if recs else 0.0:.6f}",
            "abstain_fraction": f"{(len(abstain) / len(recs)) if recs else 0.0:.6f}",
            "mean_confidence": f"{mean(confidences):.6f}",
            "mean_abs_donor_margin": f"{mean(donor_margins):.6f}",
            "decision_sources": ",".join(decision_sources),
            "positions": ",".join(map(str, positions)),
            "wrong_nonzero_positions": ",".join(map(str, wrong_nonzero)),
            "abstain_positions": ",".join(map(str, abstain)),
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "window_index",
                "n_sites",
                "n_wrong_nonzero",
                "n_abstain",
                "n_correct",
                "wrong_nonzero_fraction",
                "abstain_fraction",
                "mean_confidence",
                "mean_abs_donor_margin",
                "decision_sources",
                "positions",
                "wrong_nonzero_positions",
                "abstain_positions",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out}")
    print(f"Rows classified: {len(pos_class)}")

if __name__ == "__main__":
    main()
