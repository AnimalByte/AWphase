#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 8 ]]; then
  echo "Usage:"
  echo "  bash scripts/run_phase8_pbwt_hmm_window.sh CHROM START END SOURCE_WINDOW PANEL_BCF GENETIC_MAP TRUTH_VCF TRUTH_BED"
  exit 1
fi

CHROM="$1"
START="$2"
END="$3"
SOURCE_WINDOW="$4"
PANEL_BCF="$5"
GENETIC_MAP="$6"
TRUTH_VCF="$7"
TRUTH_BED="$8"

ROOT="results/phase7a_windows/${SOURCE_WINDOW}"
VARIANT_JSON="${ROOT}/variants.window.json"
BACKBONE="${ROOT}/local_calls.phase6c.tsv"
mkdir -p \
  "$ROOT/phase8_pbwt" \
  "$ROOT/phase8_pbwt_hmm" \
  "$ROOT/phase8_pbwt_v2" \
  "$ROOT/phase8_pbwt_hmm_v2" \
  "$ROOT/phase8_block_scaffold" \
  "$ROOT/phase8_pbwt_v3" \
  "$ROOT/phase8_pbwt_hmm_v3"
TIMING_TSV="$ROOT/phase8_timing.tsv"

printf 'label\tmethod\tstep\truntime_seconds\tstarted_epoch\tended_epoch\n' > "$TIMING_TSV"

record_timing() {
  local method="$1"
  local step="$2"
  local started="$3"
  local ended
  ended="$(date +%s)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SOURCE_WINDOW" \
    "$method" \
    "$step" \
    "$((ended - started))" \
    "$started" \
    "$ended" >> "$TIMING_TSV"
}

for f in "$PANEL_BCF" "$GENETIC_MAP" "$TRUTH_VCF" "$TRUTH_BED" "$VARIANT_JSON" "$BACKBONE"; do
  if [[ ! -s "$f" ]]; then
    echo "MISSING required file: $f" >&2
    exit 1
  fi
done

RUST_PBWT_BIN="${AWPHASE_RUST_PBWT_BIN:-target/release/phase_panel_pbwt}"
if [[ -x "$RUST_PBWT_BIN" ]]; then
  RUST_PBWT_CMD=("$RUST_PBWT_BIN")
else
  RUST_PBWT_CMD=(cargo run --release --quiet --bin phase_panel_pbwt --)
fi

RUST_SCAFFOLD_BIN="${AWPHASE_RUST_SCAFFOLD_BIN:-target/release/phase_panel_scaffold}"
if [[ -x "$RUST_SCAFFOLD_BIN" ]]; then
  RUST_SCAFFOLD_CMD=("$RUST_SCAFFOLD_BIN")
else
  RUST_SCAFFOLD_CMD=(cargo run --release --quiet --bin phase_panel_scaffold --)
fi

echo "===== ${SOURCE_WINDOW}: Phase8 PBWT projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --out-tsv "$ROOT/local_calls.phase8pbwt.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt/candidates.pbwt.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt/summary.pbwt.json" \
  > "$ROOT/phase8_pbwt/run.pbwt.log" 2>&1
record_timing "AWPhase_Phase8_PBWT" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt/summary.pbwt.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT + HMM projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt-hmm \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --genetic-map "$GENETIC_MAP" \
  --out-tsv "$ROOT/local_calls.phase8pbwt_hmm.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt_hmm/candidates.pbwt_hmm.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt_hmm/summary.pbwt_hmm.json" \
  > "$ROOT/phase8_pbwt_hmm/run.pbwt_hmm.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_HMM" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt_hmm/summary.pbwt_hmm.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_hmm.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_hmm"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT v2 exact-prefix projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt-v2 \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --out-tsv "$ROOT/local_calls.phase8pbwt_v2.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt_v2/candidates.pbwt_v2.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt_v2/summary.pbwt_v2.json" \
  > "$ROOT/phase8_pbwt_v2/run.pbwt_v2.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_V2_exact_prefix" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt_v2/summary.pbwt_v2.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_v2.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_v2"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT v2 + forward-backward HMM projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt-hmm-v2 \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --genetic-map "$GENETIC_MAP" \
  --hmm-min-margin 0.25 \
  --out-tsv "$ROOT/local_calls.phase8pbwt_hmm_v2.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt_hmm_v2/candidates.pbwt_hmm_v2.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt_hmm_v2/summary.pbwt_hmm_v2.json" \
  > "$ROOT/phase8_pbwt_hmm_v2/run.pbwt_hmm_v2.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_HMM_V2_forward_backward" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt_hmm_v2/summary.pbwt_hmm_v2.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_hmm_v2.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_hmm_v2"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT v2 + HMM block scaffolding ====="
started="$(date +%s)"
"${RUST_SCAFFOLD_CMD[@]}" \
  --panel-bcf "$PANEL_BCF" \
  --input-tsv "$ROOT/local_calls.phase8pbwt_hmm_v2.tsv" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --genetic-map "$GENETIC_MAP" \
  --out-tsv "$ROOT/local_calls.phase8pbwt_hmm_v2_scaffold.tsv" \
  --out-joins-tsv "$ROOT/phase8_block_scaffold/joins.pbwt_hmm_v2_scaffold.tsv" \
  --out-summary-json "$ROOT/phase8_block_scaffold/summary.pbwt_hmm_v2_scaffold.json" \
  --min-anchors-per-block "${AWPHASE_SCAFFOLD_MIN_ANCHORS_PER_BLOCK:-2}" \
  --max-anchors-per-block "${AWPHASE_SCAFFOLD_MAX_ANCHORS_PER_BLOCK:-48}" \
  --top-k "${AWPHASE_SCAFFOLD_TOP_K:-64}" \
  --min-donors "${AWPHASE_SCAFFOLD_MIN_DONORS:-32}" \
  --min-confidence "${AWPHASE_SCAFFOLD_MIN_CONFIDENCE:-0.55}" \
  --min-margin "${AWPHASE_SCAFFOLD_MIN_MARGIN:-8.0}" \
  --max-gap-cm "${AWPHASE_SCAFFOLD_MAX_GAP_CM:-0.50}" \
  --max-gap-bp "${AWPHASE_SCAFFOLD_MAX_GAP_BP:-2000000}" \
  > "$ROOT/phase8_block_scaffold/run.pbwt_hmm_v2_scaffold.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_HMM_V2_block_scaffold" "phase_panel_scaffold" "$started"

cat "$ROOT/phase8_block_scaffold/summary.pbwt_hmm_v2_scaffold.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_hmm_v2_scaffold.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_hmm_v2_scaffold"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT v3 bidirectional-prefix projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt-v3 \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --out-tsv "$ROOT/local_calls.phase8pbwt_v3.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt_v3/candidates.pbwt_v3.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt_v3/summary.pbwt_v3.json" \
  > "$ROOT/phase8_pbwt_v3/run.pbwt_v3.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_V3_bidirectional_prefix" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt_v3/summary.pbwt_v3.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_v3.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_v3"

echo
echo "===== ${SOURCE_WINDOW}: Phase8 PBWT v3 + forward-backward HMM projection ====="
started="$(date +%s)"
"${RUST_PBWT_CMD[@]}" \
  --mode pbwt-hmm-v3 \
  --panel-bcf "$PANEL_BCF" \
  --backbone-local-calls-tsv "$BACKBONE" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --genetic-map "$GENETIC_MAP" \
  --hmm-min-margin 0.25 \
  --out-tsv "$ROOT/local_calls.phase8pbwt_hmm_v3.tsv" \
  --out-candidates-tsv "$ROOT/phase8_pbwt_hmm_v3/candidates.pbwt_hmm_v3.tsv" \
  --out-summary-json "$ROOT/phase8_pbwt_hmm_v3/summary.pbwt_hmm_v3.json" \
  > "$ROOT/phase8_pbwt_hmm_v3/run.pbwt_hmm_v3.log" 2>&1
record_timing "AWPhase_Phase8_PBWT_HMM_V3_bidirectional_forward_backward" "phase_panel_pbwt" "$started"

cat "$ROOT/phase8_pbwt_hmm_v3/summary.pbwt_hmm_v3.json"

PYTHONPATH=python python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase8pbwt_hmm_v3.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --sample HG002 \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase8pbwt_hmm_v3"
