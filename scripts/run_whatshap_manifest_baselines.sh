#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-results/phase8f_manifests/phase8f_windows.window_beds.tsv}"
SPLIT_FILTER="${2:-all}"
REF_FASTA="${3:-references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}"
JOBS="${JOBS:-2}"

if [[ ! -s "$MANIFEST" ]]; then
  echo "MISSING manifest: $MANIFEST" >&2
  exit 1
fi

if [[ ! -s "$REF_FASTA" ]]; then
  echo "MISSING reference FASTA: $REF_FASTA" >&2
  exit 1
fi

mkdir -p results/experiments/whatshap_manifest

pacbio_bam_for_chrom() {
  local chrom="$1"
  if [[ "$chrom" == "chr6" ]]; then
    echo "data/raw/hg002_chr6/pacbio_hifi/HG002.pacbio-revio-hifi.30x.chr6_25_35mb.bam"
    return
  fi
  echo "data/raw/hg002_${chrom}/pacbio_hifi/HG002.pacbio-revio-hifi.30x.${chrom}.bam"
}

start_job() {
  while [[ "$(jobs -pr | wc -l)" -ge "$JOBS" ]]; do
    sleep 1
  done
  "$@" &
  PIDS+=("$!")
}

# shellcheck disable=SC2329
run_one() {
  local chrom="$1"
  local start="$2"
  local end="$3"
  local label="$4"
  local bam="$5"
  local variant_json="$6"
  local truth_vcf="$7"
  local truth_bed="$8"
  local method="$9"

  local out_label="${label}_whatshap_${method}"
  local metrics="baselines/whatshap/${out_label}/truth_eval.metrics.json"
  local log="results/experiments/whatshap_manifest/${out_label}.log"

  if [[ -s "$metrics" ]]; then
    echo "SKIP existing ${out_label}"
    return 0
  fi

  echo "RUN ${out_label}"
  bash scripts/run_whatshap_baseline_window.sh \
    "$chrom" \
    "$start" \
    "$end" \
    "$out_label" \
    "$bam" \
    "$variant_json" \
    "$truth_vcf" \
    "$truth_bed" \
    "$REF_FASTA" \
    > "$log" 2>&1
}

PIDS=()

while IFS=$'\t' read -r SPLIT CHROM START END LABEL _READS _SNP _BEDINT BAM _PANEL _JSON TRUTH_VCF TRUTH_BED; do
  if [[ "$SPLIT" == "split" ]]; then
    continue
  fi
  if [[ "$SPLIT_FILTER" != "all" && "$SPLIT" != "$SPLIT_FILTER" ]]; then
    continue
  fi

  ROOT="results/phase7a_windows/${LABEL}_illumina30x_phase7a"
  VARIANT_JSON="${ROOT}/variants.window.json"
  if [[ ! -s "$VARIANT_JSON" ]]; then
    echo "MISSING variant JSON for ${LABEL}: $VARIANT_JSON" >&2
    echo "Run scripts/run_phase8_pbwt_hmm_manifest.sh for this split first." >&2
    exit 1
  fi

  if [[ -s "$BAM" ]]; then
    start_job run_one "$CHROM" "$START" "$END" "$LABEL" "$BAM" "$VARIANT_JSON" "$TRUTH_VCF" "$TRUTH_BED" "illumina30x"
  else
    echo "SKIP missing Illumina BAM for ${LABEL}: $BAM"
  fi

  PACBIO_BAM="$(pacbio_bam_for_chrom "$CHROM")"
  if [[ -s "$PACBIO_BAM" ]]; then
    start_job run_one "$CHROM" "$START" "$END" "$LABEL" "$PACBIO_BAM" "$VARIANT_JSON" "$TRUTH_VCF" "$TRUTH_BED" "pacbio"
  else
    echo "SKIP missing PacBio BAM for ${LABEL}: $PACBIO_BAM"
  fi
done < "$MANIFEST"

status=0
for pid in "${PIDS[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done

exit "$status"
