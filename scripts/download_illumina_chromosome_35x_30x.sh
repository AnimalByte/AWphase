#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: bash scripts/download_illumina_chromosome_35x_30x.sh chr1"
  exit 1
fi

CHR="$1"

BAM_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"

OUTDIR="data/raw/hg002_${CHR}/illumina_chrom"
mkdir -p "$OUTDIR"

BAM35="$OUTDIR/HG002.illumina.35x.${CHR}.bam"
BAM30="$OUTDIR/HG002.illumina.30x.from35x.${CHR}.bam"

echo "===== Extract whole chromosome 35x ====="
echo "CHR=$CHR"
echo "BAM35=$BAM35"

if [[ ! -s "$BAM35" ]]; then
  samtools view -@ 8 -b "$BAM_URL" "$CHR" -o "$BAM35"
  samtools index "$BAM35"
else
  echo "Already exists: $BAM35"
fi

echo
echo "===== Downsample whole chromosome to approx 30x ====="
echo "BAM30=$BAM30"

if [[ ! -s "$BAM30" ]]; then
  samtools view -@ 8 -b -s 42.857 "$BAM35" -o "$BAM30"
  samtools index "$BAM30"
else
  echo "Already exists: $BAM30"
fi

echo
echo "===== idxstats ====="
samtools idxstats "$BAM35" | awk -v c="$CHR" '$1==c'
samtools idxstats "$BAM30" | awk -v c="$CHR" '$1==c'

echo
echo "DONE"
