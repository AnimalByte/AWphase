#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import pysam

def is_biallelic_variant(rec):
    return rec.alts is not None and len(rec.alts) == 1

def is_snp(rec):
    return is_biallelic_variant(rec) and len(rec.ref) == 1 and len(rec.alts[0]) == 1

def is_het_sample(rec, sample):
    if sample not in rec.samples:
        return False
    gt = rec.samples[sample].get("GT")
    if gt is None or len(gt) != 2 or None in gt:
        return False
    return set(gt) == {0, 1}

def count_het_variants(vcf, chrom, start, end, sample):
    total = 0
    snp = 0
    for rec in vcf.fetch(chrom, start - 1, end):
        if rec.pos < start or rec.pos > end:
            continue
        if not is_biallelic_variant(rec):
            continue
        if not is_het_sample(rec, sample):
            continue
        total += 1
        if is_snp(rec):
            snp += 1
    return total, snp

def n_fraction(ref, chrom, start, end):
    if ref is None:
        return ""
    try:
        seq = ref.fetch(chrom, start - 1, end).upper()
    except Exception:
        return ""
    if not seq:
        return 1.0
    return seq.count("N") / len(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--reference-fasta", default="")
    ap.add_argument("--window-bp", type=int, default=5_000_000)
    ap.add_argument("--step-bp", type=int, default=5_000_000)
    ap.add_argument("--edge-buffer-bp", type=int, default=10_000_000)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-snp-variants", type=int, default=500)
    ap.add_argument("--min-reads", type=int, default=50_000)
    ap.add_argument("--max-n-frac", type=float, default=0.01)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    if not bam.has_index():
        raise SystemExit(f"BAM is not indexed: {args.bam}")

    vcf = pysam.VariantFile(args.truth_vcf)
    if args.chrom not in vcf.header.contigs:
        raise SystemExit(f"{args.chrom} not found in VCF header contigs.")

    ref = None
    if args.reference_fasta and Path(args.reference_fasta).exists():
        ref = pysam.FastaFile(args.reference_fasta)

    chrom_len = bam.get_reference_length(args.chrom)

    rows = []
    start_min = max(1, args.edge_buffer_bp)
    end_max = chrom_len - args.edge_buffer_bp

    start = start_min
    while start + args.window_bp - 1 <= end_max:
        end = start + args.window_bp - 1

        reads = bam.count(args.chrom, start - 1, end)
        total_vars, snp_vars = count_het_variants(vcf, args.chrom, start, end, args.sample)
        nf = n_fraction(ref, args.chrom, start, end) if ref else ""

        pass_basic = (
            reads >= args.min_reads
            and snp_vars >= args.min_snp_variants
            and (nf == "" or float(nf) <= args.max_n_frac)
        )

        score = reads * max(1, snp_vars)

        rows.append({
            "chrom": args.chrom,
            "start": start,
            "end": end,
            "window_bp": args.window_bp,
            "reads": reads,
            "het_variant_count": total_vars,
            "snp_variant_count": snp_vars,
            "n_fraction": nf,
            "pass_basic": int(pass_basic),
            "score": score,
        })

        start += args.step_bp

    bam.close()
    vcf.close()
    if ref:
        ref.close()

    rows.sort(key=lambda r: (int(r["pass_basic"]), int(r["snp_variant_count"]), int(r["reads"])), reverse=True)

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = ["chrom", "start", "end", "window_bp", "reads", "het_variant_count", "snp_variant_count", "n_fraction", "pass_basic", "score"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows[:args.top_n]:
            w.writerow(r)

    print({
        "windows_scored": len(rows),
        "written": min(len(rows), args.top_n),
        "out": str(out),
    })

if __name__ == "__main__":
    main()
