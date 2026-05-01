#!/usr/bin/env python3
import csv
from pathlib import Path

in_manifest = Path("results/phase8f_manifests/phase8f_windows.tsv")
out_manifest = Path("results/phase8f_manifests/phase8f_windows.window_beds.tsv")

rows = []

with open(in_manifest) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        chrom = r["chrom"]
        start = int(r["start"])
        end = int(r["end"])
        label = r["label"]

        full_bed = Path(r["truth_bed"])
        if not full_bed.exists():
            r["split"] = "exclude_missing_truth_bed"
            rows.append(r)
            continue

        # BED is 0-based half-open. Your runner window is 1-based inclusive.
        q0 = start - 1
        q1 = end

        # Existing convention:
        # HG002_GRCh38_v5.0q_smvar.chr1.benchmark.bed
        # -> HG002_GRCh38_v5.0q_smvar.chr1_20_25mb.benchmark.bed
        name = full_bed.name
        old = f".{chrom}.benchmark.bed"
        if name.endswith(old):
            out_name = name.replace(old, f".{label}.benchmark.bed")
        else:
            out_name = f"{full_bed.stem}.{label}.bed"

        out_bed = full_bed.parent / out_name

        n = 0
        with open(full_bed) as inp, open(out_bed, "w") as out:
            for line in inp:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                c, s, e = parts[0], int(parts[1]), int(parts[2])
                if c != chrom:
                    continue
                if e <= q0 or s >= q1:
                    continue

                parts[1] = str(max(s, q0))
                parts[2] = str(min(e, q1))
                if int(parts[2]) > int(parts[1]):
                    out.write("\t".join(parts) + "\n")
                    n += 1

        r["truth_bed"] = str(out_bed)
        r["benchmark_bed_intervals"] = str(n)

        if r["split"] == "train" and n == 0:
            r["split"] = "exclude_no_window_truth_bed"

        rows.append(r)

fields = [
    "split", "chrom", "start", "end", "label",
    "reads", "snp_variant_count", "benchmark_bed_intervals",
    "bam", "panel", "json", "truth_vcf", "truth_bed",
]

with open(out_manifest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print("wrote", out_manifest)
for r in rows:
    print(r["split"], r["label"], "bed_intervals=", r["benchmark_bed_intervals"], "bed=", r["truth_bed"])
