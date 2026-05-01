#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def load_tsv(path):
    rows = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for row in reader:
            pos = None
            for key in ("pos", "site_pos", "position", "variant_pos"):
                if key in row and str(row[key]).strip():
                    pos = int(str(row[key]).strip())
                    break
            if pos is not None:
                rows[pos] = row
    return rows, fieldnames

def pick_value(row, preferred_names, contains_all=None):
    for k in preferred_names:
        if k in row and str(row[k]).strip():
            return str(row[k]).strip(), k

    if contains_all:
        for k, v in row.items():
            lk = k.lower()
            if all(token in lk for token in contains_all) and str(v).strip():
                return str(v).strip(), k

    return "", ""

def get_truth(row):
    return pick_value(
        row,
        [
            "truth_phase_state",
            "phase_state_truth",
            "truth_state",
            "truth_phase",
            "truth",
        ],
        contains_all=["truth", "phase"],
    )

def get_pred_local(row):
    return pick_value(
        row,
        [
            "local_phase_state",
            "pred_phase_state",
            "predicted_phase_state",
            "phase_state_pred",
            "pred_state",
        ],
        contains_all=["pred", "phase"],
    )

def get_pred_stitched(row):
    return pick_value(
        row,
        [
            "stitched_phase_state",
            "pred_phase_state",
            "predicted_phase_state",
            "phase_state_pred",
            "pred_state",
        ],
        contains_all=["pred", "phase"],
    )

def classify(pred, truth):
    if not pred or not truth:
        return "unknown"
    if pred == "0":
        return "abstain"
    if pred == truth:
        return "correct"
    return "wrong"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", required=True)
    ap.add_argument("--stitched", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    local, local_fields = load_tsv(args.local)
    stitched, stitched_fields = load_tsv(args.stitched)

    print("LOCAL_FIELDS:", local_fields)
    print("STITCHED_FIELDS:", stitched_fields)

    rows = []
    for pos in sorted(set(local) & set(stitched)):
        l = local[pos]
        s = stitched[pos]

        truth, truth_col_l = get_truth(l)
        truth2, truth_col_s = get_truth(s)
        if not truth:
            truth = truth2

        lpred, lpred_col = get_pred_local(l)
        spred, spred_col = get_pred_stitched(s)

        lclass = classify(lpred, truth)
        sclass = classify(spred, truth)

        pred_changed = (lpred != spred)

        if lclass == "unknown" or sclass == "unknown":
            delta = "unknown"
        elif lclass == "correct" and sclass == "wrong":
            delta = "bridge_hurt"
        elif lclass == "wrong" and sclass == "correct":
            delta = "bridge_helped"
        elif lclass == "correct" and sclass == "abstain":
            delta = "bridge_lost_call"
        elif lclass == "abstain" and sclass == "correct":
            delta = "bridge_recovered_call"
        elif pred_changed:
            delta = "pred_changed_same_correctness"
        else:
            delta = "same_bucket"

        rows.append({
            "pos": pos,
            "truth_phase_state": truth,
            "local_pred": lpred,
            "local_pred_col": lpred_col,
            "local_class": lclass,
            "stitched_pred": spred,
            "stitched_pred_col": spred_col,
            "stitched_class": sclass,
            "pred_changed": pred_changed,
            "delta": delta,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "pos",
                "truth_phase_state",
                "local_pred",
                "local_pred_col",
                "local_class",
                "stitched_pred",
                "stitched_pred_col",
                "stitched_class",
                "pred_changed",
                "delta",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out}")

if __name__ == "__main__":
    main()
