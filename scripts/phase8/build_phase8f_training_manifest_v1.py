#!/usr/bin/env python3
import csv
from pathlib import Path

# Reserve these for later external holdout. Do not train on them.
HOLDOUT = {
    ("chr1", 100000000, 104999999),
    ("chr1", 165000000, 169999999),
    ("chr22", 25000000, 29999999),
    ("chr22", 30000000, 34999999),
    ("chr22", 35000000, 39999999),
}

# Already-used six core windows; keep them in training for Phase8F.
CORE = {
    ("chr1", 15000000, 19999999),
    ("chr20", 25000000, 29999999),
    ("chr20", 30000000, 34999999),
    ("chr20", 45000000, 49999999),
    ("chr22", 15000000, 19999999),
    ("chr22", 20000000, 24999999),
}

window_files = [
    ("chr1",  Path("results/window_selection/chr1.illumina30x.fullchrom.top_windows.with_bed.tsv")),
    ("chr20", Path("results/window_selection/chr20.fulljson.top_windows.tsv")),
    ("chr22", Path("results/window_selection/chr22.illumina30x.fullchrom.top_windows.with_bed.tsv")),
]

# Known resource paths.
resources = {
    "chr1": {
        "bam": "data/raw/hg002_chr1/illumina_chrom/HG002.illumina.30x.from35x.chr1.bam",
        "panel": "data/panels/hgdp_1kg_hg38/full/hgdp1kgp_chr1.filtered.SNV_INDEL.phased.shapeit5.bcf",
        "json": "data/derived/full_chrom_variants/hg002_chr1_full.v5q.variants.real.json",
        "truth_vcf": "data/truth/hg002_chr1/HG002_GRCh38_v5.0q_smvar.chr1.vcf.gz",
        "truth_bed": "data/truth/hg002_chr1/HG002_GRCh38_v5.0q_smvar.chr1.benchmark.bed",
    },
    "chr20": {
        "bam": "baselines/illumina_30x/HG002.illumina.chr20.30x.bam",
        "panel": "data/panels/hgdp_1kg_hg38/full/hgdp1kgp_chr20.filtered.SNV_INDEL.phased.shapeit5.bcf",
        "json": "data/derived/full_chrom_variants/hg002_chr20_full.v5q.variants.real.json",
        "truth_vcf": "data/truth/hg002_chr20/HG002_GRCh38_v5.0q_smvar.chr20.vcf.gz",
        "truth_bed": "data/truth/hg002_chr20/HG002_GRCh38_v5.0q_smvar.chr20.benchmark.bed",
    },
    "chr22": {
        "bam": "data/raw/hg002_chr22/illumina_chrom/HG002.illumina.30x.from35x.chr22.bam",
        "panel": "data/panels/hgdp_1kg_hg38/full/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf",
        "json": "data/derived/full_chrom_variants/hg002_chr22_full.v5q.variants.real.json",
        "truth_vcf": "data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz",
        "truth_bed": "data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.benchmark.bed",
    },
}

rows = []

for chrom, path in window_files:
    if not path.exists():
        print("MISSING", path)
        continue

    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                start = int(r["start"])
                end = int(r["end"])
                reads = int(float(r.get("reads", 0)))
                snps = int(float(r.get("snp_variant_count", r.get("variant_count", 0))))
            except Exception:
                continue

            key = (chrom, start, end)
            label = f"{chrom}_{start//1000000}_{(end+1)//1000000}mb"

            if key in HOLDOUT:
                split = "holdout_reserved"
            else:
                split = "train"

            # Keep only usable training windows.
            bed_intervals = int(float(r.get("benchmark_bed_intervals", 1)))
            if split == "train":
                if reads < 50000 or snps < 500 or bed_intervals == 0:
                    continue

            rr = {
                "split": split,
                "chrom": chrom,
                "start": start,
                "end": end,
                "label": label,
                "reads": reads,
                "snp_variant_count": snps,
                "benchmark_bed_intervals": bed_intervals,
            }
            rr.update(resources[chrom])
            rows.append(rr)

# Deduplicate by chrom/start/end. Prefer core windows and train over holdout.
dedup = {}
for r in rows:
    key = (r["chrom"], r["start"], r["end"])
    old = dedup.get(key)
    if old is None:
        dedup[key] = r
    else:
        if (r["chrom"], r["start"], r["end"]) in CORE:
            dedup[key] = r

rows = list(dedup.values())

# Sort: train first, then holdout, then chrom/start.
rows.sort(key=lambda r: (r["split"] != "train", r["chrom"], int(r["start"])))

out = Path("results/phase8f_manifests/phase8f_windows.tsv")
out.parent.mkdir(parents=True, exist_ok=True)

fields = [
    "split", "chrom", "start", "end", "label",
    "reads", "snp_variant_count", "benchmark_bed_intervals",
    "bam", "panel", "json", "truth_vcf", "truth_bed",
]

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print("wrote", out)
print("train windows:", sum(1 for r in rows if r["split"] == "train"))
print("holdout reserved:", sum(1 for r in rows if r["split"] == "holdout_reserved"))
