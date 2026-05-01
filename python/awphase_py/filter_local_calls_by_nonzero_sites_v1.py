#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
    return 0

def read_nonzero_positions(path):
    pos = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                p = int(float(r["pos"]))
            except Exception:
                continue
            if state(r) != 0:
                pos.add(p)
    return pos

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--reference-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--mode", choices=["keep_reference_nonzero", "keep_intersection_nonzero"], default="keep_reference_nonzero")
    args = ap.parse_args()

    ref_pos = read_nonzero_positions(args.reference_tsv)
    input_nonzero = read_nonzero_positions(args.input_tsv)

    if args.mode == "keep_intersection_nonzero":
        keep = ref_pos & input_nonzero
    else:
        keep = ref_pos

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    kept = 0
    nonzero_kept = 0

    with open(args.input_tsv) as infh, open(out, "w", newline="") as outfh:
        reader = csv.DictReader(infh, delimiter="\t")
        fields = reader.fieldnames
        writer = csv.DictWriter(outfh, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for r in reader:
            n += 1
            p = int(float(r["pos"]))
            if p not in keep:
                # Preserve row but force abstain, so evaluator sees same row universe.
                if "local_phase_state" in r:
                    r["local_phase_state"] = 0
                if "phase_state" in r:
                    r["phase_state"] = 0
                if "block_id" in r:
                    r["block_id"] = "unassigned"
            else:
                kept += 1
                if state(r) != 0:
                    nonzero_kept += 1
            writer.writerow(r)

    print({
        "input": args.input_tsv,
        "reference": args.reference_tsv,
        "out": args.out_tsv,
        "mode": args.mode,
        "rows": n,
        "reference_nonzero_positions": len(ref_pos),
        "input_nonzero_positions": len(input_nonzero),
        "keep_positions": len(keep),
        "rows_in_keep_positions": kept,
        "nonzero_rows_kept": nonzero_kept,
    })

if __name__ == "__main__":
    main()
