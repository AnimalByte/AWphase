#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample-col", type=int, default=9)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with opener(args.vcf) as fh, open(out, "w") as outfh:
        outfh.write("pos\tblock_id\tphase_state\tlocal_phase_state\tconfidence\n")
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos = fields[1]
            fmt = fields[8].split(":")
            sample = fields[args.sample_col].split(":")
            fmt_map = {k: i for i, k in enumerate(fmt)}

            gt = sample[fmt_map["GT"]] if "GT" in fmt_map and fmt_map["GT"] < len(sample) else "./."
            ps = sample[fmt_map["PS"]] if "PS" in fmt_map and fmt_map["PS"] < len(sample) else "unassigned"

            if gt == "0|1":
                phase_state = "1"
            elif gt == "1|0":
                phase_state = "-1"
            else:
                phase_state = "0"

            block_id = f"PS_{ps}" if ps not in (".", "", "unassigned") else "unassigned"
            outfh.write(f"{pos}\t{block_id}\t{phase_state}\t{phase_state}\t1.0\n")

if __name__ == "__main__":
    main()
