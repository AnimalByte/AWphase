#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import pysam


def find_meta_cols(fieldnames):
    sample_keys = ["sample", "sample_id", "s", "id", "donor_id"]
    group_keys = ["group", "superpopulation", "super_pop", "population", "pop", "continental_pop", "genetic_region"]
    sample_col = next((k for k in sample_keys if k in fieldnames), None)
    group_col = next((k for k in group_keys if k in fieldnames), None)
    return sample_col, group_col


def load_group_map(meta_tsv):
    if not meta_tsv:
        return {}
    with open(meta_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sample_col, group_col = find_meta_cols(reader.fieldnames)
        if sample_col is None or group_col is None:
            return {}
        out = {}
        for row in reader:
            s = str(row[sample_col]).strip()
            g = str(row[group_col]).strip()
            if s:
                out[s] = g if g else "UNK"
        return out


def gt_is_het_01(gt):
    if gt is None or len(gt) != 2:
        return False
    a, b = gt
    if a is None or b is None:
        return False
    return {a, b} == {0, 1}


def allele_code(a):
    if a is None:
        return 0
    if a == 0:
        return -1
    if a == 1:
        return 1
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-vcf", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--region", required=True)
    ap.add_argument("--meta-tsv", default="")
    ap.add_argument("--out-variants-json", required=True)
    ap.add_argument("--out-site-priors-on-json", required=True)
    ap.add_argument("--out-site-priors-off-json", required=True)
    ap.add_argument("--out-panel-candidates-json", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    chrom, rest = args.region.split(":")
    start, end = [int(x.replace(",", "")) for x in rest.split("-")]

    target_vcf = pysam.VariantFile(args.target_vcf)
    panel_bcf = pysam.VariantFile(args.panel_bcf)
    group_map = load_group_map(args.meta_tsv)

    matched_sites = []
    site_priors = []
    hap_arrays = defaultdict(list)

    target_seen = 0
    kept = 0
    skipped = defaultdict(int)

    panel_samples = list(panel_bcf.header.samples)

    for rec in target_vcf.fetch(chrom, start - 1, end):
        target_seen += 1

        if len(rec.alts or []) != 1:
            skipped["target_non_biallelic"] += 1
            continue

        ref = rec.ref
        alt = rec.alts[0]

        if len(ref) != 1 or len(alt) != 1:
            skipped["target_non_snv"] += 1
            continue

        if args.sample not in rec.samples:
            skipped["target_missing_sample"] += 1
            continue

        gt = rec.samples[args.sample].get("GT")
        if not gt_is_het_01(gt):
            skipped["target_not_het01"] += 1
            continue

        panel_match = None
        for prec in panel_bcf.fetch(chrom, rec.pos - 1, rec.pos):
            if prec.pos == rec.pos and prec.ref == ref and len(prec.alts or []) == 1 and prec.alts[0] == alt:
                panel_match = prec
                break

        if panel_match is None:
            skipped["no_exact_panel_match"] += 1
            continue

        kept += 1
        matched_sites.append({
            "pos": rec.pos,
            "ref_allele": ref,
            "alt_allele": alt,
            "genotype": [0, 1],
        })

        alt_count = 0
        total_count = 0

        for s in panel_samples:
            gt = panel_match.samples[s].get("GT")
            if gt is None or len(gt) != 2:
                a0 = a1 = 0
            else:
                a0 = allele_code(gt[0])
                a1 = allele_code(gt[1])

                for a in gt:
                    if a is not None and a in (0, 1):
                        total_count += 1
                        if a == 1:
                            alt_count += 1

            hap_arrays[(s, 0)].append(a0)
            hap_arrays[(s, 1)].append(a1)

        af = (alt_count / total_count) if total_count > 0 else 0.5
        bias = 2.0 * af - 1.0
        site_priors.append({"pos": rec.pos, "bias": bias})

    candidates = []
    for sample in panel_samples:
        group = group_map.get(sample, "UNK")
        candidates.append({
            "donor_id": sample,
            "hap_id": 0,
            "group": group,
            "alleles": hap_arrays[(sample, 0)],
        })
        candidates.append({
            "donor_id": sample,
            "hap_id": 1,
            "group": group,
            "alleles": hap_arrays[(sample, 1)],
        })

    Path(args.out_variants_json).parent.mkdir(parents=True, exist_ok=True)
    for path, obj in [
        (args.out_variants_json, matched_sites),
        (args.out_panel_candidates_json, candidates),
        (args.out_site_priors_on_json, {"entries": site_priors}),
        (args.out_site_priors_off_json, {"entries": site_priors}),
    ]:
        with open(path, "w") as fh:
            json.dump(obj, fh)

    summary = {
        "region": args.region,
        "target_records_seen": target_seen,
        "matched_sites_kept": kept,
        "skip_counts": dict(skipped),
        "panel_samples": len(panel_samples),
        "candidate_haplotypes": len(candidates),
    }
    with open(args.out_summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
