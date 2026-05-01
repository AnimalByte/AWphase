#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def fval(x, default=0.0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return float(x)
    except Exception:
        return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    rows = []
    with open(args.audit_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            correctness = r.get("correctness", "")
            if correctness not in {"correct", "wrong"}:
                continue

            support = fval(r.get("panel_support"))
            conflict = fval(r.get("panel_conflict"))
            margin = fval(r.get("panel_margin"))
            conf = fval(r.get("panel_confidence"))
            second = fval(r.get("best_vs_second_margin"))
            samples = fval(r.get("panel_samples"))
            haps = fval(r.get("panel_haplotypes"))
            anchors = fval(r.get("anchors"))
            n_blocks = fval(r.get("n_block_candidates"))
            orient_sites = fval(r.get("orientation_sites"))
            orient_margin = fval(r.get("orientation_margin"))

            rows.append({
                "label": 1 if correctness == "correct" else 0,
                "pos": r.get("pos", ""),
                "panel_confidence": conf,
                "panel_margin": margin,
                "panel_support": support,
                "panel_conflict": conflict,
                "panel_log_support_conflict_ratio": (support + 1.0) / (conflict + 1.0),
                "panel_samples": samples,
                "panel_haplotypes": haps,
                "anchors": anchors,
                "best_vs_second_margin": second,
                "n_block_candidates": n_blocks,
                "orientation_sites": orient_sites,
                "orientation_margin": orient_margin,
            })

    fields = [
        "label", "pos",
        "panel_confidence", "panel_margin", "panel_support", "panel_conflict",
        "panel_log_support_conflict_ratio",
        "panel_samples", "panel_haplotypes", "anchors",
        "best_vs_second_margin", "n_block_candidates",
        "orientation_sites", "orientation_margin",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print({"rows": len(rows), "positives": sum(int(r["label"]) for r in rows)})

if __name__ == "__main__":
    main()
