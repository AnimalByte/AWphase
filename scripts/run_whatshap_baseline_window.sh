#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 9 ]]; then
  echo "Usage:"
  echo "  bash scripts/run_whatshap_baseline_window.sh CHROM START END LABEL BAM VARIANT_JSON TRUTH_VCF TRUTH_BED REF_FASTA"
  exit 1
fi

CHROM="$1"
START="$2"
END="$3"
LABEL="$4"
BAM="$5"
VARIANT_JSON="$6"
TRUTH_VCF="$7"
TRUTH_BED="$8"
REF_FASTA="$9"

SAMPLE="HG002"
OUT="baselines/whatshap/${LABEL}"
mkdir -p "$OUT"
TIMING_TSV="$OUT/run_timing.tsv"

case "$LABEL" in
  *_whatshap_illumina30x) METHOD="WhatsHap_Illumina30x" ;;
  *_whatshap_pacbio) METHOD="WhatsHap_PacBio_HiFi30x" ;;
  *) METHOD="WhatsHap" ;;
esac

printf 'label\tmethod\tstep\truntime_seconds\tstarted_epoch\tended_epoch\n' > "$TIMING_TSV"

record_timing() {
  local method="$1"
  local step="$2"
  local started="$3"
  local ended
  ended="$(date +%s)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LABEL" \
    "$method" \
    "$step" \
    "$((ended - started))" \
    "$started" \
    "$ended" >> "$TIMING_TSV"
}

UNPHASED="$OUT/unphased.input.vcf.gz"
RAW="$OUT/whatshap.raw.vcf"
PHASED="$OUT/whatshap.phased.vcf.gz"

echo "===== WhatsHap baseline inputs ====="
echo "CHROM=$CHROM"
echo "START=$START"
echo "END=$END"
echo "LABEL=$LABEL"
echo "BAM=$BAM"
echo "VARIANT_JSON=$VARIANT_JSON"
echo "TRUTH_VCF=$TRUTH_VCF"
echo "TRUTH_BED=$TRUTH_BED"
echo "REF_FASTA=$REF_FASTA"

for f in "$BAM" "$VARIANT_JSON" "$TRUTH_VCF" "$TRUTH_BED" "$REF_FASTA"; do
  if [[ ! -s "$f" ]]; then
    echo "MISSING required file: $f"
    exit 1
  fi
done

if [[ ! -s "${REF_FASTA}.fai" ]]; then
  samtools faidx "$REF_FASTA"
fi

if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
  samtools index "$BAM"
fi

echo
echo "===== Build unphased VCF input ====="
bcftools view -r "${CHROM}:${START}-${END}" -s "$SAMPLE" "$TRUTH_VCF" \
  | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {gsub(/\|/,"/",$10); print}' \
  | bgzip -c > "$UNPHASED"

tabix -f -p vcf "$UNPHASED"

echo
echo "===== Run WhatsHap ====="
started="$(date +%s)"
whatshap phase \
  --reference "$REF_FASTA" \
  --sample "$SAMPLE" \
  --ignore-read-groups \
  --output "$RAW" \
  "$UNPHASED" \
  "$BAM" \
  > "$OUT/whatshap.stdout.log" \
  2> "$OUT/whatshap.stderr.log"
record_timing "$METHOD" "whatshap_phase" "$started"

echo
echo "===== Compress phased VCF ====="
bcftools view -Oz -o "$PHASED" "$RAW"
tabix -f -p vcf "$PHASED"

echo
echo "===== Convert WhatsHap VCF to local_calls.tsv ====="
PYTHONPATH=python python python/awphase_py/whatshap_vcf_to_local_calls_v1.py \
  --phased-vcf "$PHASED" \
  --variant-json "$VARIANT_JSON" \
  --out-tsv "$OUT/local_calls.tsv" \
  --sample "$SAMPLE" \
  > "$OUT/convert.log" 2>&1

cat "$OUT/convert.log"

echo
echo "===== Evaluate WhatsHap baseline ====="
python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$OUT/local_calls.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "$OUT/truth_eval"

echo
echo "===== Metrics ====="
cat "$OUT/truth_eval.metrics.json"
