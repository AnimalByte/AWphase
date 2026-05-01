#!/usr/bin/env python3
import argparse
import pysam
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--sample", default="HG002")
    args = ap.parse_args()

    inv = pysam.VariantFile(args.in_vcf)
    sample = args.sample if args.sample in inv.header.samples else list(inv.header.samples)[0]

    out_path = Path(args.out_vcf)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    out = pysam.VariantFile(str(out_path), "w", header=inv.header)

    n = 0
    n_het = 0
    n_unphased = 0

    for rec in inv.fetch():
        if sample not in rec.samples:
            continue

        sd = rec.samples[sample]
        gt = sd.get("GT")

        if gt is not None and len(gt) == 2 and None not in gt:
            if tuple(gt) in {(0, 1), (1, 0)}:
                n_het += 1
            try:
                sd.phased = False
                n_unphased += 1
            except Exception:
                pass

        # Remove PS if present so WhatsHap creates fresh phase sets.
        if "PS" in sd:
            try:
                del sd["PS"]
            except Exception:
                pass

        out.write(rec)
        n += 1

    out.close()
    inv.close()

    print({
        "input": args.in_vcf,
        "output": args.out_vcf,
        "sample": sample,
        "records_written": n,
        "het_records": n_het,
        "records_unphased": n_unphased,
    })

if __name__ == "__main__":
    main()
