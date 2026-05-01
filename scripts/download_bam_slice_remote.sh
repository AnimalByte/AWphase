#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 5 ]]; then
  echo "Usage:"
  echo "  bash scripts/download_bam_slice_remote.sh BAM_URL REGION OUT_BAM THREADS SAMPLE_FRAC"
  echo
  echo "Examples:"
  echo "  SAMPLE_FRAC=1.0 for no downsampling"
  echo "  SAMPLE_FRAC=0.857 for 35x -> approx 30x"
  exit 1
fi

BAM_URL="$1"
REGION="$2"
OUT_BAM="$3"
THREADS="$4"
SAMPLE_FRAC="$5"

mkdir -p "$(dirname "$OUT_BAM")"

echo "BAM_URL=$BAM_URL"
echo "REGION=$REGION"
echo "OUT_BAM=$OUT_BAM"
echo "SAMPLE_FRAC=$SAMPLE_FRAC"

if [[ "$SAMPLE_FRAC" == "1.0" || "$SAMPLE_FRAC" == "1" ]]; then
  samtools view -@ "$THREADS" -b "$BAM_URL" "$REGION" \
    -o "$OUT_BAM"
else
  # Fixed seed 42 + fraction.
  samtools view -@ "$THREADS" -b -s "42${SAMPLE_FRAC}" "$BAM_URL" "$REGION" \
    -o "$OUT_BAM"
fi

samtools index "$OUT_BAM"

echo
echo "===== idxstats ====="
CHROM="${REGION%%:*}"
samtools idxstats "$OUT_BAM" | awk -v c="$CHROM" '$1==c || NR<=3'

echo
echo "DONE: $OUT_BAM"
