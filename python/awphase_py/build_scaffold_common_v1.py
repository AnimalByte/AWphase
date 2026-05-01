#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pysam


def load_variant_json(path):
    with open(path) as fh:
        obj = json.load(fh)
    if isinstance(obj, dict) and "variants" in obj and isinstance(obj["variants"], list):
        return obj, obj["variants"]
    if isinstance(obj, list):
        return None, obj
    raise ValueError("Unsupported variant JSON structure")


def dump_variant_json(wrapper, variants, out_path):
    with open(out_path, "w") as fh:
        if wrapper is None:
            json.dump(variants, fh)
        else:
            wrapper = dict(wrapper)
            wrapper["variants"] = variants
            json.dump(wrapper, fh)


def load_reads_json(path):
    with open(path) as fh:
        obj = json.load(fh)
    if not isinstance(obj, list):
        raise ValueError("Expected reads JSON to be a list of observations")
    return obj


def get_af_from_record(rec):
    if "AF" in rec.info:
        af = rec.info["AF"]
        if isinstance(af, (list, tuple)):
            return float(af[0])
        return float(af)

    # Fallback: compute AF from GTs
    alt_alleles = 0
    total_alleles = 0
    for s in rec.samples.values():
        gt = s.get("GT")
        if gt is None:
            continue
        for a in gt:
            if a is None:
                continue
            total_alleles += 1
            if a == 1:
                alt_alleles += 1
    if total_alleles == 0:
        return None
    return alt_alleles / total_alleles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-variants-json", required=True)
    ap.add_argument("--target-reads-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("--min-af", type=float, default=0.05)
    ap.add_argument("--max-af", type=float, default=0.95)
    ap.add_argument("--out-variants-json", required=True)
    ap.add_argument("--out-reads-json", required=True)
    ap.add_argument("--out-positions-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    wrapper, variants = load_variant_json(args.target_variants_json)
    target_positions = {int(v["pos"]) for v in variants if isinstance(v, dict) and "pos" in v}

    panel = pysam.VariantFile(args.panel_bcf)

    scaffold_positions = set()
    panel_records_seen = 0
    panel_records_matching_target = 0
    non_biallelic = 0
    non_snv = 0
    af_missing = 0
    af_outside = 0

    for rec in panel.fetch(args.chrom):
        panel_records_seen += 1

        if rec.pos not in target_positions:
            continue
        panel_records_matching_target += 1

        if len(rec.alts or []) != 1:
            non_biallelic += 1
            continue

        ref = rec.ref
        alt = rec.alts[0]
        if len(ref) != 1 or len(alt) != 1:
            non_snv += 1
            continue

        af = get_af_from_record(rec)
        if af is None:
            af_missing += 1
            continue

        if af < args.min_af or af > args.max_af:
            af_outside += 1
            continue

        scaffold_positions.add(int(rec.pos))

    filtered_variants = [
        v for v in variants
        if isinstance(v, dict) and "pos" in v and int(v["pos"]) in scaffold_positions
    ]

    reads = load_reads_json(args.target_reads_json)
    filtered_reads = [
        r for r in reads
        if isinstance(r, dict) and "site_pos" in r and int(r["site_pos"]) in scaffold_positions
    ]

    Path(args.out_variants_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_reads_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_positions_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)

    dump_variant_json(wrapper, filtered_variants, args.out_variants_json)

    with open(args.out_reads_json, "w") as fh:
        json.dump(filtered_reads, fh)

    with open(args.out_positions_tsv, "w") as fh:
        for pos in sorted(scaffold_positions):
            fh.write(f"{args.chrom}\t{pos}\n")

    summary = {
        "chrom": args.chrom,
        "target_variant_positions": len(target_positions),
        "panel_records_seen": panel_records_seen,
        "panel_records_matching_target": panel_records_matching_target,
        "common_scaffold_positions": len(scaffold_positions),
        "filtered_variants": len(filtered_variants),
        "filtered_read_observations": len(filtered_reads),
        "min_af": args.min_af,
        "max_af": args.max_af,
        "skip_counts": {
            "non_biallelic": non_biallelic,
            "non_snv": non_snv,
            "af_missing": af_missing,
            "af_outside": af_outside,
        },
    }

    with open(args.out_summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
