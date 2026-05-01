#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict, Counter

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

def is_filled(row):
    return str(row.get("phase7a_panel_filled", "0")).strip().lower() in {"1", "true", "t", "yes"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--fill-candidates-tsv", required=True)
    ap.add_argument("--site-comparison-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-tsv", required=True)
    ap.add_argument("--min-orientation-sites", type=int, default=2)
    args = ap.parse_args()

    local_rows = read_tsv(args.local_calls_tsv)
    fill_rows = read_tsv(args.fill_candidates_tsv)
    comp_rows = read_tsv(args.site_comparison_tsv)

    local_by_pos = {}
    for r in local_rows:
        pos = ival(r.get("pos"), None)
        if pos is not None:
            local_by_pos[pos] = r

    fill_by_pos = {}
    for r in fill_rows:
        pos = ival(r.get("pos"), None)
        if pos is not None:
            fill_by_pos[pos] = r

    comp_by_pos = {}
    for r in comp_rows:
        pos = ival(r.get("pos"), None)
        if pos is not None:
            comp_by_pos[pos] = r

    # Learn block orientation using non-panel-filled sites only.
    # This avoids circularly using the filled site to orient itself.
    orient_counts = defaultdict(lambda: {"same": 0, "opposite": 0})

    for pos, cr in comp_by_pos.items():
        lr = local_by_pos.get(pos, {})
        if is_filled(lr):
            continue

        block = cr.get("block_id", lr.get("block_id", ""))
        pred = ival(cr.get("pred_state"), 0)
        truth = ival(cr.get("truth_state"), 0)
        matched = str(cr.get("matched_exact_truth", "1")).strip()

        if not block or block.lower() in {"unassigned", "none", "na", "."}:
            continue
        if pred == 0 or truth == 0:
            continue
        if matched in {"0", "false", "False"}:
            continue

        if pred == truth:
            orient_counts[block]["same"] += 1
        elif pred == -truth:
            orient_counts[block]["opposite"] += 1

    block_orientation = {}
    for block, c in orient_counts.items():
        same = c["same"]
        opp = c["opposite"]
        n = same + opp

        if n < args.min_orientation_sites:
            block_orientation[block] = {
                "orientation": "",
                "orientation_sites": n,
                "same": same,
                "opposite": opp,
                "orientation_margin": abs(same - opp),
                "orientation_reason": "too_few_sites",
            }
        else:
            o = 1 if same >= opp else -1
            block_orientation[block] = {
                "orientation": o,
                "orientation_sites": n,
                "same": same,
                "opposite": opp,
                "orientation_margin": abs(same - opp),
                "orientation_reason": "nonfilled_truth_sites",
            }

    out_rows = []
    counts = Counter()

    for pos, lr in sorted(local_by_pos.items()):
        if not is_filled(lr):
            continue

        counts["filled_sites"] += 1

        fr = fill_by_pos.get(pos, {})
        cr = comp_by_pos.get(pos, {})

        block = cr.get("block_id", lr.get("block_id", ""))
        pred = ival(cr.get("pred_state", lr.get("local_phase_state", lr.get("phase_state"))), 0)
        truth = ival(cr.get("truth_state"), 0)
        matched = str(cr.get("matched_exact_truth", "1")).strip()

        bo = block_orientation.get(block, {
            "orientation": "",
            "orientation_sites": 0,
            "same": 0,
            "opposite": 0,
            "orientation_margin": 0,
            "orientation_reason": "no_orientation_sites",
        })

        orientation = bo["orientation"]
        oriented_pred = ""

        if pred == 0 or truth == 0 or matched in {"0", "false", "False"}:
            correctness = "no_truth_comparison"
            counts["filled_no_truth_comparison"] += 1
        elif orientation == "":
            correctness = "unknown_no_block_orientation"
            counts["filled_unknown_no_block_orientation"] += 1
        else:
            oriented_pred = pred * int(orientation)
            if oriented_pred == truth:
                correctness = "correct"
                counts["filled_correct_oriented"] += 1
            else:
                correctness = "wrong"
                counts["filled_wrong_oriented"] += 1

        row = {
            "pos": pos,
            "block_id": block,
            "raw_pred_state": pred,
            "truth_state": truth,
            "block_orientation": orientation,
            "oriented_pred_state": oriented_pred,
            "correctness": correctness,
            "orientation_sites": bo["orientation_sites"],
            "orientation_same": bo["same"],
            "orientation_opposite": bo["opposite"],
            "orientation_margin": bo["orientation_margin"],
            "orientation_reason": bo["orientation_reason"],
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
        out_rows.append(row)

    fields = [
        "pos", "block_id", "raw_pred_state", "truth_state",
        "block_orientation", "oriented_pred_state", "correctness",
        "orientation_sites", "orientation_same", "orientation_opposite",
        "orientation_margin", "orientation_reason",
        "panel_confidence", "panel_margin", "panel_support", "panel_conflict",
        "panel_samples", "panel_haplotypes", "anchors",
        "best_vs_second_margin", "n_block_candidates",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    summary_rows = []
    for k in [
        "filled_sites",
        "filled_correct_oriented",
        "filled_wrong_oriented",
        "filled_no_truth_comparison",
        "filled_unknown_no_block_orientation",
    ]:
        summary_rows.append({"metric": k, "value": counts.get(k, 0)})

    wrong = counts.get("filled_wrong_oriented", 0)
    correct = counts.get("filled_correct_oriented", 0)
    denom = correct + wrong

    summary_rows.append({
        "metric": "oriented_fill_accuracy",
        "value": (correct / denom if denom else ""),
    })

    with open(args.out_summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"], delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    print({r["metric"]: r["value"] for r in summary_rows})

if __name__ == "__main__":
    main()
