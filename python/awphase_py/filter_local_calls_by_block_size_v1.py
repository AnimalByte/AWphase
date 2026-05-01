#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from collections import defaultdict

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
    return 0

def set_state(row, s):
    if "local_phase_state" in row:
        row["local_phase_state"] = int(s)
    if "phase_state" in row:
        row["phase_state"] = int(s)

def block_id(row):
    b = str(row.get("block_id", "")).strip()
    if b == "" or b.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return b

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-block-sites", type=int, required=True)
    ap.add_argument("--count-nonzero-only", action="store_true", default=True)
    args = ap.parse_args()

    rows = []
    blocks = defaultdict(list)

    with open(args.input_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        for r in reader:
            rows.append(r)
            s = state(r)
            b = block_id(r)
            if b and s != 0:
                blocks[b].append(r)

    block_sizes = {b: len(v) for b, v in blocks.items()}

    kept_nonzero = 0
    masked_nonzero = 0

    out_rows = []
    for r in rows:
        rr = dict(r)
        s = state(rr)
        b = block_id(rr)

        if s != 0:
            if not b or block_sizes.get(b, 0) < args.min_block_sites:
                set_state(rr, 0)
                rr["block_id"] = "unassigned"
                masked_nonzero += 1
            else:
                kept_nonzero += 1

        out_rows.append(rr)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    print({
        "input": args.input_tsv,
        "out": args.out_tsv,
        "min_block_sites": args.min_block_sites,
        "blocks_total": len(blocks),
        "blocks_kept": sum(1 for _, n in block_sizes.items() if n >= args.min_block_sites),
        "kept_nonzero": kept_nonzero,
        "masked_nonzero": masked_nonzero,
    })

if __name__ == "__main__":
    main()
