#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import pysam

def is_het_01(gt):
    if gt is None or len(gt) != 2 or None in gt:
        return False
    return set(gt) == {0, 1}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--snps-only", action="store_true")
    args = ap.parse_args()

    vf = pysam.VariantFile(args.truth_vcf)

    if args.sample not in vf.header.samples:
        raise SystemExit(f"Sample {args.sample} not found. Samples={list(vf.header.samples)[:10]}...")

    out = []
    skipped_multiallelic = 0
    skipped_nonhet = 0
    skipped_nonsnp = 0

    for rec in vf.fetch(args.chrom, args.start - 1, args.end):
        if rec.pos < args.start or rec.pos > args.end:
            continue
        if rec.alts is None or len(rec.alts) != 1:
            skipped_multiallelic += 1
            continue

        gt = rec.samples[args.sample].get("GT")
        if not is_het_01(gt):
            skipped_nonhet += 1
            continue

        ref = str(rec.ref).upper()
        alt = str(rec.alts[0]).upper()

        if args.snps_only and not (len(ref) == 1 and len(alt) == 1):
            skipped_nonsnp += 1
            continue

        # Important: genotype is unphased/sorted to avoid truth-phase leakage.
        out.append({
            "chrom": args.chrom,
            "pos": int(rec.pos),
            "ref_allele": ref,
            "alt_allele": alt,
            "genotype": [0, 1],
        })

    vf.close()

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2) + "\n")

    print({
        "truth_vcf": args.truth_vcf,
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "out_json": args.out_json,
        "rows": len(out),
        "skipped_multiallelic": skipped_multiallelic,
        "skipped_nonhet": skipped_nonhet,
        "skipped_nonsnp": skipped_nonsnp,
        "snps_only": args.snps_only,
    })

if __name__ == "__main__":
    main()
