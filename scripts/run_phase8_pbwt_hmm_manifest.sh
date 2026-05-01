#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-results/phase8f_manifests/phase8f_windows.window_beds.tsv}"
SPLIT_FILTER="${2:-train}"
READSET="${AWPHASE_READSET:-illumina30x}"

case "$READSET" in
  illumina30x)
    SOURCE_SUFFIX="${AWPHASE_SOURCE_SUFFIX:-_illumina30x_phase7a}"
    SUMMARY_OUT_DIR="${AWPHASE_PHASE8_SUMMARY_DIR:-results/phase8/pbwt_hmm_manifest}"
    PANEL_RELIANCE_OUT_DIR="${AWPHASE_PANEL_RELIANCE_DIR:-results/phase8/panel_reliance_manifest}"
    CALIBRATION_OUT_DIR="${AWPHASE_CALIBRATION_DIR:-results/phase8/confidence_calibration_manifest}"
    ;;
  pacbio | hifi | pacbio_hifi30x)
    READSET="pacbio_hifi30x"
    SOURCE_SUFFIX="${AWPHASE_SOURCE_SUFFIX:-_pacbio_hifi30x_mcs128_phase7a}"
    SUMMARY_OUT_DIR="${AWPHASE_PHASE8_SUMMARY_DIR:-results/phase8/pbwt_hmm_pacbio_hifi_mcs128_manifest}"
    PANEL_RELIANCE_OUT_DIR="${AWPHASE_PANEL_RELIANCE_DIR:-results/phase8/panel_reliance_pacbio_hifi_mcs128_manifest}"
    CALIBRATION_OUT_DIR="${AWPHASE_CALIBRATION_DIR:-results/phase8/confidence_calibration_pacbio_hifi_mcs128_manifest}"
    export AWPHASE_MAX_COMPONENT_SITES="${AWPHASE_MAX_COMPONENT_SITES:-128}"
    export AWPHASE_LOCAL_REFINE_ITERS="${AWPHASE_LOCAL_REFINE_ITERS:-20}"
    ;;
  *)
    SOURCE_SUFFIX="${AWPHASE_SOURCE_SUFFIX:-_${READSET}_phase7a}"
    SUMMARY_OUT_DIR="${AWPHASE_PHASE8_SUMMARY_DIR:-results/phase8/pbwt_hmm_${READSET}_manifest}"
    PANEL_RELIANCE_OUT_DIR="${AWPHASE_PANEL_RELIANCE_DIR:-results/phase8/panel_reliance_${READSET}_manifest}"
    CALIBRATION_OUT_DIR="${AWPHASE_CALIBRATION_DIR:-results/phase8/confidence_calibration_${READSET}_manifest}"
    ;;
esac

if [[ ! -s "$MANIFEST" ]]; then
  echo "MISSING manifest: $MANIFEST" >&2
  exit 1
fi

mkdir -p results/experiments/phase8_pbwt_hmm_manifest

pacbio_bam_for_chrom() {
  local chrom="$1"
  local start="$2"
  local end="$3"

  if [[ "$chrom" == "chr6" && "$start" -ge 25000000 && "$end" -le 35000000 ]]; then
    echo "data/raw/hg002_chr6/pacbio_hifi/HG002.pacbio-revio-hifi.30x.chr6_25_35mb.bam"
    return
  fi

  echo "data/raw/hg002_${chrom}/pacbio_hifi/HG002.pacbio-revio-hifi.30x.${chrom}.bam"
}

bam_for_readset() {
  local chrom="$1"
  local start="$2"
  local end="$3"
  local manifest_bam="$4"

  case "$READSET" in
    illumina30x)
      echo "$manifest_bam"
      ;;
    pacbio_hifi30x)
      pacbio_bam_for_chrom "$chrom" "$start" "$end"
      ;;
    *)
      if [[ -n "${AWPHASE_BAM_OVERRIDE:-}" ]]; then
        echo "$AWPHASE_BAM_OVERRIDE"
      else
        echo "unsupported AWPHASE_READSET=$READSET; set AWPHASE_BAM_OVERRIDE" >&2
        return 1
      fi
      ;;
  esac
}

tail -n +2 "$MANIFEST" |
while IFS=$'\t' read -r SPLIT CHROM START END LABEL _READS _SNP _BEDINT BAM PANEL JSON TRUTH_VCF TRUTH_BED; do
  if [[ "$SPLIT_FILTER" != "all" && "$SPLIT" != "$SPLIT_FILTER" ]]; then
    continue
  fi

  SOURCE_WINDOW="${LABEL}${SOURCE_SUFFIX}"
  ROOT="results/phase7a_windows/${SOURCE_WINDOW}"
  MAP="data/maps/beagle_grch38/no_chr_in_chrom_field/plink.${CHROM}.GRCh38.map"
  INPUT_BAM="$(bam_for_readset "$CHROM" "$START" "$END" "$BAM")"

  if [[ ! -s "$INPUT_BAM" ]]; then
    echo "MISSING ${READSET} BAM for ${SOURCE_WINDOW}: $INPUT_BAM" >&2
    exit 1
  fi

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
      "$INPUT_BAM" \
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
  --split "$SPLIT_FILTER" \
  --out-dir "$SUMMARY_OUT_DIR" \
  --source-suffix "$SOURCE_SUFFIX"

PYTHONPATH=python python scripts/phase8/analyze_panel_reliance_v1.py \
  --manifest "$MANIFEST" \
  --split "$SPLIT_FILTER" \
  --out-dir "$PANEL_RELIANCE_OUT_DIR" \
  --source-suffix "$SOURCE_SUFFIX"

PYTHONPATH=python python scripts/phase8/calibrate_phase8_confidence_v1.py \
  --manifest "$MANIFEST" \
  --split "$SPLIT_FILTER" \
  --out-dir "$CALIBRATION_OUT_DIR" \
  --source-suffix "$SOURCE_SUFFIX"
