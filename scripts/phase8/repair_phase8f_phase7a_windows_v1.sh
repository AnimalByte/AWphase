#!/usr/bin/env bash
set +e +u
set +o pipefail

MANIFEST="results/phase8f_manifests/phase8f_windows.window_beds.tsv"
mkdir -p results/experiments

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r SPLIT CHR START END LABEL READS SNP BEDINT BAM PANEL JSON TRUTH_VCF TRUTH_BED; do
  [[ "$SPLIT" == "train" ]] || continue

  
  # Strip possible Windows carriage returns from manifest fields.
  SPLIT="${SPLIT//$'\r'/}"
  CHR="${CHR//$'\r'/}"
  START="${START//$'\r'/}"
  END="${END//$'\r'/}"
  LABEL="${LABEL//$'\r'/}"
  READS="${READS//$'\r'/}"
  SNP="${SNP//$'\r'/}"
  BEDINT="${BEDINT//$'\r'/}"
  BAM="${BAM//$'\r'/}"
  PANEL="${PANEL//$'\r'/}"
  JSON="${JSON//$'\r'/}"
  TRUTH_VCF="${TRUTH_VCF//$'\r'/}"
  TRUTH_BED="${TRUTH_BED//$'\r'/}"


  ROOT="results/phase7a_windows/${LABEL}_illumina30x_phase7a"

  if [[ -s "$ROOT/local_calls.phase7a.tsv" && -s "$ROOT/fill_candidates.phase7a.tsv" && -s "$ROOT/truth_eval_phase7a.site_comparison.oriented.tsv" ]]; then
    echo "SKIP complete Phase7A: $LABEL"
    continue
  fi

  echo
  echo "================================================================="
  echo "REPAIR Phase7A: $LABEL"
  echo "CHR=$CHR START=$START END=$END"
  echo "TRUTH_BED=$TRUTH_BED"
  echo "================================================================="

  bash scripts/run_phase7a_window.sh \
    "$CHR" \
    "$START" \
    "$END" \
    "${LABEL}_illumina30x_phase7a" \
    "$BAM" \
    "$PANEL" \
    "$JSON" \
    "$TRUTH_VCF" \
    "$TRUTH_BED" \
    > "results/experiments/phase8f_repair_phase7a_${LABEL}.log" 2>&1

  rc=$?
  echo "$LABEL rc=$rc"
  tail -n 80 "results/experiments/phase8f_repair_phase7a_${LABEL}.log"

  if [[ "$rc" != "0" ]]; then
    echo "FAILED $LABEL; continuing"
  fi
done
