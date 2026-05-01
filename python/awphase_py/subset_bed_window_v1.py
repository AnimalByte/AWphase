#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-bed", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)  # 1-based inclusive
    ap.add_argument("--end", type=int, required=True)    # 1-based inclusive
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    # BED is 0-based half-open. Convert requested interval.
    win0 = args.start - 1
    win1 = args.end

    out = Path(args.out_bed)
    out.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(args.in_bed) as fh, open(out, "w") as oh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            s = int(parts[1])
            e = int(parts[2])
            if chrom != args.chrom:
                continue
            os = max(s, win0)
            oe = min(e, win1)
            if os < oe:
                parts[1] = str(os)
                parts[2] = str(oe)
                oh.write("\t".join(parts) + "\n")
                n += 1

    print({"out_bed": str(out), "intervals": n})

if __name__ == "__main__":
    main()
