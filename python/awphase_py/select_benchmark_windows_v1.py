#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
import pysam

def load_variants(path, chrom):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    out = []
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue

        c = str(r.get("chrom", r.get("contig", chrom)))
        if c and c != chrom:
            continue

        ref = str(r.get("ref_allele", r.get("ref", "")))
        alt = str(r.get("alt_allele", r.get("alt", "")))
        vtype = "snp" if len(ref) == 1 and len(alt) == 1 else "indel_or_complex"

        out.append((pos, vtype))

    out.sort()
    return out

def count_variants_in_window(variants, start, end):
    snp = 0
    total = 0
    for pos, typ in variants:
        if pos < start:
            continue
        if pos > end:
            break
        total += 1
        if typ == "snp":
            snp += 1
    return total, snp

def n_fraction(ref, chrom, start, end):
    if ref is None:
        return ""
    seq = ref.fetch(chrom, start - 1, end).upper()
    if not seq:
        return 1.0
    return seq.count("N") / len(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--reference-fasta", default="")
    ap.add_argument("--window-bp", type=int, default=5_000_000)
    ap.add_argument("--step-bp", type=int, default=5_000_000)
    ap.add_argument("--edge-buffer-bp", type=int, default=5_000_000)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-snp-variants", type=int, default=500)
    ap.add_argument("--min-reads", type=int, default=50_000)
    ap.add_argument("--max-n-frac", type=float, default=0.01)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    if not bam.has_index():
        raise SystemExit(f"BAM is not indexed: {args.bam}")

    chrom_len = bam.get_reference_length(args.chrom)
    variants = load_variants(args.variant_json, args.chrom)

    ref = None
    if args.reference_fasta and Path(args.reference_fasta).exists():
        ref = pysam.FastaFile(args.reference_fasta)

    rows = []

    start_min = max(1, args.edge_buffer_bp)
    end_max = chrom_len - args.edge_buffer_bp

    start = start_min
    while start + args.window_bp - 1 <= end_max:
        end = start + args.window_bp - 1

        total_vars, snp_vars = count_variants_in_window(variants, start, end)

        try:
            reads = bam.count(args.chrom, start - 1, end)
        except Exception:
            reads = 0

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
            "variant_count": total_vars,
            "snp_variant_count": snp_vars,
            "n_fraction": nf,
            "pass_basic": int(pass_basic),
            "score": score,
        })

        start += args.step_bp

    bam.close()
    if ref:
        ref.close()

    rows.sort(key=lambda r: (int(r["pass_basic"]), int(r["score"])), reverse=True)

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = ["chrom", "start", "end", "window_bp", "reads", "variant_count", "snp_variant_count", "n_fraction", "pass_basic", "score"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows[:args.top_n]:
            w.writerow(r)

    print({"windows_scored": len(rows), "written": min(len(rows), args.top_n), "out": str(out)})

if __name__ == "__main__":
    main()
