#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

# Some AWPhase TSV rows can carry large annotation/audit fields.
# Raise Python's default CSV parser field limit.
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fill-candidates-tsv", required=True)
    ap.add_argument("--oriented-site-comparison-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    comp = read_tsv(args.oriented_site_comparison_tsv)
    fills = read_tsv(args.fill_candidates_tsv)

    truth_by_pos = {}
    orient_by_block = {}

    for r in comp:
        pos = ival(r.get("pos"), None)
        if pos is not None:
            truth_by_pos[pos] = ival(r.get("truth_state"), 0)

        block = str(r.get("block_id", "")).strip()
        orient = str(r.get("block_orientation", "")).strip()
        if block and orient not in {"", "0", "None"}:
            orient_by_block[block] = ival(orient, 0)

    rows = []

    for r in fills:
        pos = ival(r.get("pos"), None)
        if pos is None:
            continue

        block = str(r.get("block_id", "")).strip()
        pred = ival(r.get("pred_state"), 0)
        truth = truth_by_pos.get(pos, 0)
        orient = orient_by_block.get(block, 0)

        if pred == 0 or truth == 0 or orient == 0:
            continue

        oriented_pred = pred * orient
        label = 1 if oriented_pred == truth else 0

        support = fval(r.get("support"))
        conflict = fval(r.get("conflict"))
        margin = fval(r.get("margin"))
        confidence = fval(r.get("confidence"))

        rows.append({
            "label": label,
            "source_window": args.source_window,
            "pos": pos,
            "block_id": block,
            "candidate_accepted_by_rule": ival(r.get("accepted"), 0),
            "pred_state": pred,
            "truth_state": truth,
            "block_orientation": orient,
            "oriented_pred_state": oriented_pred,
            "panel_confidence": confidence,
            "panel_margin": margin,
            "panel_support": support,
            "panel_conflict": conflict,
            "support_conflict_ratio": (support + 1.0) / (conflict + 1.0),
            "panel_samples": fval(r.get("panel_samples")),
            "panel_haplotypes": fval(r.get("panel_haplotypes")),
            "anchors": fval(r.get("anchors")),
            "best_vs_second_margin": fval(r.get("best_vs_second_margin")),
            "n_block_candidates": fval(r.get("n_block_candidates")),
        })

    fields = [
        "label", "source_window", "pos", "block_id",
        "candidate_accepted_by_rule",
        "pred_state", "truth_state", "block_orientation", "oriented_pred_state",
        "panel_confidence", "panel_margin", "panel_support", "panel_conflict",
        "support_conflict_ratio", "panel_samples", "panel_haplotypes",
        "anchors", "best_vs_second_margin", "n_block_candidates",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    positives = sum(int(r["label"]) for r in rows)
    print({
        "source_window": args.source_window,
        "rows": len(rows),
        "positives": positives,
        "negatives": len(rows) - positives,
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
