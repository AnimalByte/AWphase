#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def read_positions(path):
    pos = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                pos.add(int(float(r["pos"])))
            except Exception:
                pass
    return pos

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--calls-tsv", action="append", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--chrom-start", type=int, default=None)
    ap.add_argument("--chrom-end", type=int, default=None)
    args = ap.parse_args()

    keep_pos = set()
    for p in args.calls_tsv:
        keep_pos |= read_positions(p)

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_out = 0

    with open(args.obs_tsv) as infh, open(out, "w", newline="") as outfh:
        reader = csv.DictReader(infh, delimiter="\t")
        writer = csv.DictWriter(outfh, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for r in reader:
            n_in += 1
            try:
                p = int(float(r.get("site_pos", r.get("pos", 0))))
            except Exception:
                continue

            if args.chrom_start is not None and p < args.chrom_start:
                continue
            if args.chrom_end is not None and p > args.chrom_end:
                continue
            if p not in keep_pos:
                continue

            writer.writerow(r)
            n_out += 1

    print({
        "obs_in": n_in,
        "obs_out": n_out,
        "positions_kept": len(keep_pos),
        "out_tsv": str(out),
    })

if __name__ == "__main__":
    main()
