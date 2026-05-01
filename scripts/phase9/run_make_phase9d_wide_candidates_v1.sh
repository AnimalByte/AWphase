#!/usr/bin/env bash
set +e +u
set +o pipefail

MANIFEST="results/phase8f_manifests/phase8f_windows.window_beds.tsv"
mkdir -p results/phase9/logs

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r SPLIT CHR START END LABEL _READS _SNP _BEDINT _BAM PANEL _JSON _TRUTH_VCF _TRUTH_BED; do
  SPLIT="${SPLIT//$'\r'/}"
  CHR="${CHR//$'\r'/}"
  START="${START//$'\r'/}"
  END="${END//$'\r'/}"
  LABEL="${LABEL//$'\r'/}"
  PANEL="${PANEL//$'\r'/}"

  [[ "$SPLIT" == "train" ]] || continue

  ROOT="results/phase7a_windows/${LABEL}_illumina30x_phase7a"
  GMAP="data/maps/beagle_grch38/no_chr_in_chrom_field/plink.${CHR}.GRCh38.map"
  OUT="$ROOT/phase9d_block_bridge_candidates.wide.tsv"
  LOG="results/phase9/logs/make_phase9d_wide_${LABEL}.log"

  echo
  echo "===== Phase9D wide bridge candidates: $LABEL ====="

  missing=0
  for f in \
    "$ROOT/local_calls.phase6c.tsv" \
    "$ROOT/truth_eval_phase6c.site_comparison.tsv" \
    "$ROOT/variants.window.json" \
    "$PANEL" \
    "$GMAP"
  do
    if [[ ! -s "$f" ]]; then
      echo "MISSING prerequisite: $f"
      missing=1
    fi
  done

  if [[ "$missing" == "1" ]]; then
    echo "SKIP $LABEL"
    continue
  fi

  PYTHONPATH=python python python/awphase_py/phase9/make_phase9a_block_bridge_candidates_v1.py \
    --source-window "${LABEL}_illumina30x_phase7a" \
    --local-calls-tsv "$ROOT/local_calls.phase6c.tsv" \
    --truth-comparison-tsv "$ROOT/truth_eval_phase6c.site_comparison.tsv" \
    --variant-json "$ROOT/variants.window.json" \
    --panel-bcf "$PANEL" \
    --genetic-map "$GMAP" \
    --chrom "$CHR" \
    --start "$START" \
    --end "$END" \
    --max-gap-bp 500000 \
    --max-next-blocks 20 \
    --max-anchors-each-block 16 \
    --min-anchors-each-block 2 \
    --min-orientation-sites 2 \
    --decay-cm 0.05 \
    --top-k 64 \
    --out-tsv "$OUT" \
    > "$LOG" 2>&1

  rc=$?
  echo "$LABEL rc=$rc"
  tail -n 40 "$LOG"

  if [[ -s "$OUT" ]]; then
    echo -n "candidate rows: "
    tail -n +2 "$OUT" | wc -l
  else
    echo "NO OUTPUT: $OUT"
  fi
done
