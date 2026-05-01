#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-results/phase8f_manifests/phase8f_windows.window_beds.tsv}"
SPLIT_FILTER="${2:-train}"

if [[ ! -s "$MANIFEST" ]]; then
  echo "MISSING manifest: $MANIFEST" >&2
  exit 1
fi

mkdir -p results/experiments/phase8_pbwt_hmm_manifest

tail -n +2 "$MANIFEST" |
while IFS=$'\t' read -r SPLIT CHROM START END LABEL _READS _SNP _BEDINT BAM PANEL JSON TRUTH_VCF TRUTH_BED; do
  if [[ "$SPLIT_FILTER" != "all" && "$SPLIT" != "$SPLIT_FILTER" ]]; then
    continue
  fi

  SOURCE_WINDOW="${LABEL}_illumina30x_phase7a"
  ROOT="results/phase7a_windows/${SOURCE_WINDOW}"
  MAP="data/maps/beagle_grch38/no_chr_in_chrom_field/plink.${CHROM}.GRCh38.map"

  echo
  echo "===== ${SOURCE_WINDOW}: ensure Phase6C/7A backbone ====="
  if [[ ! -s "$ROOT/local_calls.phase6c.tsv" \
    || ! -s "$ROOT/local_calls.phase7a.tsv" \
    || ! -s "$ROOT/truth_eval_phase7a.metrics.json" \
    || ! -s "$ROOT/variants.window.json" ]]; then
    bash scripts/run_phase7a_window.sh \
      "$CHROM" \
      "$START" \
      "$END" \
      "$SOURCE_WINDOW" \
      "$BAM" \
      "$PANEL" \
      "$JSON" \
      "$TRUTH_VCF" \
      "$TRUTH_BED" \
      > "results/experiments/phase8_pbwt_hmm_manifest/phase7a_${SOURCE_WINDOW}.log" 2>&1
  else
    echo "existing backbone found: $ROOT"
  fi

  echo "===== ${SOURCE_WINDOW}: run PBWT/PBWT+HMM ====="
  bash scripts/run_phase8_pbwt_hmm_window.sh \
    "$CHROM" \
    "$START" \
    "$END" \
    "$SOURCE_WINDOW" \
    "$PANEL" \
    "$MAP" \
    "$TRUTH_VCF" \
    "$TRUTH_BED" \
    > "results/experiments/phase8_pbwt_hmm_manifest/phase8_${SOURCE_WINDOW}.log" 2>&1

  tail -n 50 "results/experiments/phase8_pbwt_hmm_manifest/phase8_${SOURCE_WINDOW}.log"
done

PYTHONPATH=python python scripts/phase8/summarize_phase8_pbwt_hmm_v1.py \
  --manifest "$MANIFEST" \
  --split "$SPLIT_FILTER"
