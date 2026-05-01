#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def pick_col(fields, candidates):
    lower = {f.lower(): f for f in fields}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None

def truthy(x):
    s = str(x).strip().lower()
    return s in {"1", "true", "t", "yes", "y", "correct"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--fill-candidates-tsv", required=True)
    ap.add_argument("--site-comparison-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-tsv", required=True)
    args = ap.parse_args()

    local = read_tsv(args.local_calls_tsv)
    fills = read_tsv(args.fill_candidates_tsv)
    comp = read_tsv(args.site_comparison_tsv)

    local_by_pos = {}
    for r in local:
        try:
            pos = int(float(r["pos"]))
        except Exception:
            continue
        local_by_pos[pos] = r

    fills_by_pos = {}
    for r in fills:
        try:
            pos = int(float(r["pos"]))
        except Exception:
            continue
        fills_by_pos[pos] = r

    comp_fields = comp[0].keys() if comp else []
    pos_col = pick_col(comp_fields, ["pos", "position"])
    correct_col = pick_col(comp_fields, [
        "correct",
        "is_correct",
        "matches_truth",
        "truth_match",
        "phase_correct",
        "is_match",
    ])
    pred_col = pick_col(comp_fields, [
        "pred_state",
        "local_phase_state",
        "phase_state",
        "pred_phase_state",
    ])
    truth_col = pick_col(comp_fields, [
        "truth_state",
        "truth_phase_state",
        "truth_local_phase_state",
    ])

    comp_by_pos = {}
    for r in comp:
        try:
            pos = int(float(r[pos_col]))
        except Exception:
            continue
        comp_by_pos[pos] = r

    rows = []
    accepted = 0
    filled = 0
    filled_with_comp = 0
    filled_correct = 0
    filled_wrong = 0

    for pos, lr in sorted(local_by_pos.items()):
        if str(lr.get("phase7a_panel_filled", "0")).strip() not in {"1", "true", "True"}:
            continue

        filled += 1
        fr = fills_by_pos.get(pos, {})
        cr = comp_by_pos.get(pos, {})

        is_correct = ""

        if cr:
            filled_with_comp += 1
            if correct_col:
                is_correct = int(truthy(cr.get(correct_col, "")))
            elif pred_col and truth_col:
                try:
                    is_correct = int(int(float(cr.get(pred_col, 0))) == int(float(cr.get(truth_col, 0))))
                except Exception:
                    is_correct = ""

        if is_correct == 1:
            filled_correct += 1
        elif is_correct == 0:
            filled_wrong += 1

        row = {
            "pos": pos,
            "phase7a_state": lr.get("local_phase_state", lr.get("phase_state", "")),
            "block_id": lr.get("block_id", ""),
            "is_correct": is_correct,
            "panel_confidence": fr.get("confidence", lr.get("phase7a_panel_confidence", "")),
            "panel_margin": fr.get("margin", lr.get("phase7a_panel_margin", "")),
            "panel_support": fr.get("support", ""),
            "panel_conflict": fr.get("conflict", ""),
            "panel_samples": fr.get("panel_samples", lr.get("phase7a_panel_samples", "")),
            "panel_haplotypes": fr.get("panel_haplotypes", lr.get("phase7a_panel_haplotypes", "")),
            "anchors": fr.get("anchors", lr.get("phase7a_panel_anchors", "")),
            "best_vs_second_margin": fr.get("best_vs_second_margin", ""),
            "n_block_candidates": fr.get("n_block_candidates", ""),
        }
        rows.append(row)

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "pos", "phase7a_state", "block_id", "is_correct",
        "panel_confidence", "panel_margin", "panel_support", "panel_conflict",
        "panel_samples", "panel_haplotypes", "anchors",
        "best_vs_second_margin", "n_block_candidates",
    ]

    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = [
        {"metric": "filled_sites", "value": filled},
        {"metric": "filled_sites_with_site_comparison", "value": filled_with_comp},
        {"metric": "filled_correct_detected", "value": filled_correct},
        {"metric": "filled_wrong_detected", "value": filled_wrong},
        {"metric": "correct_column_detected", "value": correct_col or ""},
        {"metric": "pred_column_detected", "value": pred_col or ""},
        {"metric": "truth_column_detected", "value": truth_col or ""},
    ]

    with open(args.out_summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"], delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    print({r["metric"]: r["value"] for r in summary})

if __name__ == "__main__":
    main()
