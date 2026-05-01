#!/usr/bin/env python3
import argparse
import json
from collections import defaultdict
from pathlib import Path

import pysam


def normalize_chrom_name(bam, chrom: str) -> str:
    refs = set(bam.references)
    if chrom in refs:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]
    if ("chr" + chrom) in refs:
        return "chr" + chrom
    raise ValueError(f"Chromosome {chrom!r} not found in BAM header")


def gt_is_diploid_het(gt):
    if gt is None or len(gt) != 2:
        return False
    a, b = gt
    if a is None or b is None:
        return False
    return {a, b} == {0, 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--region", required=True, help="e.g. chr20:1-5000000")
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    chrom, coords = args.region.split(":")
    start_s, end_s = coords.split("-")
    region_start = int(start_s.replace(",", ""))
    region_end = int(end_s.replace(",", ""))

    bam = pysam.AlignmentFile(args.bam, "rb")
    vcf = pysam.VariantFile(args.vcf)
    fasta = pysam.FastaFile(args.fasta)

    bam_chrom = normalize_chrom_name(bam, chrom)

    observations = []
    informative_reads = set()

    total_records = 0
    usable_het_biallelic_snvs = 0
    skipped_non_het = 0
    skipped_non_biallelic = 0
    skipped_non_snv = 0
    skipped_ref_mismatch = 0
    skipped_no_sample = 0

    obs_per_site = defaultdict(int)

    for rec in vcf.fetch(chrom, region_start - 1, region_end):
        total_records += 1

        if len(rec.alts or []) != 1:
            skipped_non_biallelic += 1
            continue

        ref = rec.ref.upper()
        alt = rec.alts[0].upper()

        if len(ref) != 1 or len(alt) != 1:
            skipped_non_snv += 1
            continue

        if args.sample not in rec.samples:
            skipped_no_sample += 1
            continue

        gt = rec.samples[args.sample].get("GT")
        if not gt_is_diploid_het(gt):
            skipped_non_het += 1
            continue

        fasta_ref = fasta.fetch(chrom, rec.pos - 1, rec.pos).upper()
        if fasta_ref != ref:
            skipped_ref_mismatch += 1
            continue

        usable_het_biallelic_snvs += 1

        for col in bam.pileup(
            bam_chrom,
            rec.pos - 1,
            rec.pos,
            truncate=True,
            min_base_quality=0,
            stepper="samtools",
            fastafile=fasta,
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            if col.reference_pos != rec.pos - 1:
                continue

            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue

                aln = pr.alignment
                if aln.is_unmapped or aln.is_duplicate or aln.is_secondary or aln.is_supplementary:
                    continue

                if aln.mapping_quality < args.min_mapq:
                    continue

                qpos = pr.query_position
                if qpos is None:
                    continue

                base = aln.query_sequence[qpos].upper()
                baseq = int(aln.query_qualities[qpos]) if aln.query_qualities is not None else 0

                if baseq < args.min_baseq:
                    continue

                if base == ref:
                    allele = -1
                elif base == alt:
                    allele = 1
                else:
                    continue

                observations.append(
                    {
                        "read_id": aln.query_name,
                        "site_pos": rec.pos,
                        "allele": allele,
                        "baseq": baseq,
                        "mapq": int(aln.mapping_quality),
                    }
                )
                informative_reads.add(aln.query_name)
                obs_per_site[rec.pos] += 1

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as fh:
        json.dump(observations, fh)

    summary = {
        "region": args.region,
        "sample": args.sample,
        "bam": args.bam,
        "vcf": args.vcf,
        "fasta": args.fasta,
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "total_vcf_records_seen": total_records,
        "usable_het_biallelic_snvs": usable_het_biallelic_snvs,
        "skipped_non_het": skipped_non_het,
        "skipped_non_biallelic": skipped_non_biallelic,
        "skipped_non_snv": skipped_non_snv,
        "skipped_ref_mismatch": skipped_ref_mismatch,
        "skipped_no_sample": skipped_no_sample,
        "observations_emitted": len(observations),
        "informative_reads": len(informative_reads),
        "sites_with_observations": len(obs_per_site),
        "mean_observations_per_observed_site": (
            sum(obs_per_site.values()) / len(obs_per_site) if obs_per_site else 0.0
        ),
    }

    out_summary = Path(args.out_summary)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    with open(out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
