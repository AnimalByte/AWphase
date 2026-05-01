#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 9 ]]; then
  echo "Usage:"
  echo "  bash scripts/run_phase7a_window.sh CHROM START END LABEL BAM PANEL_BCF VARIANT_JSON TRUTH_VCF TRUTH_BED"
  exit 1
fi

CHROM="$1"
START="$2"
END="$3"
LABEL="$4"
BAM="$5"
PANEL="$6"
VARIANT_JSON_FULL="$7"
TRUTH_VCF="$8"
TRUTH_BED="$9"

ROOT="results/phase7a_windows/${LABEL}"
mkdir -p "$ROOT"

VARIANT_JSON="$ROOT/variants.window.json"
TIMING_TSV="$ROOT/run_timing.tsv"
RUN_STARTED_EPOCH="$(date +%s)"

MIN_MAPQ="${AWPHASE_MIN_MAPQ:-20}"
MIN_BASEQ="${AWPHASE_MIN_BASEQ:-15}"
MIN_SITES_PER_FRAGMENT="${AWPHASE_MIN_SITES_PER_FRAGMENT:-2}"
MAX_EXACT_SITES="${AWPHASE_MAX_EXACT_SITES:-18}"
MAX_COMPONENT_SITES="${AWPHASE_MAX_COMPONENT_SITES:-1024}"
LOCAL_REFINE_ITERS="${AWPHASE_LOCAL_REFINE_ITERS:-8}"

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

echo "===== ${LABEL}: subset variants ====="
PYTHONPATH=python python python/awphase_py/subset_variant_json_window_v1.py \
  --in-json "$VARIANT_JSON_FULL" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --out-json "$VARIANT_JSON"

echo
echo "===== ${LABEL}: Phase6C BAM-template fragments ====="
PYTHONPATH=python python python/awphase_py/build_template_fragments_from_bam_v1.py \
  --bam "$BAM" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --min-mapq "$MIN_MAPQ" \
  --min-baseq "$MIN_BASEQ" \
  --min-sites-per-fragment "$MIN_SITES_PER_FRAGMENT" \
  --out-tsv "$ROOT/fragments.tsv" \
  --out-summary-json "$ROOT/fragments.summary.json" \
  > "$ROOT/build_fragments.log" 2>&1

tail -n 40 "$ROOT/build_fragments.log"

echo
echo "===== ${LABEL}: Phase6C solve ====="
SOLVER_ARGS=(
  --fragments-tsv "$ROOT/fragments.tsv" \
  --variant-json "$VARIANT_JSON" \
  --max-exact-sites "$MAX_EXACT_SITES" \
  --max-component-sites "$MAX_COMPONENT_SITES" \
  --local-refine-iters "$LOCAL_REFINE_ITERS" \
  --out-local-calls-tsv "$ROOT/local_calls.phase6c.tsv" \
  --out-components-tsv "$ROOT/components.phase6c.tsv" \
  --out-summary-json "$ROOT/solve.phase6c.summary.json"
)

: > "$ROOT/solve_phase6c.log"

RUST_SOLVER_BIN="${AWPHASE_RUST_SOLVER_BIN:-target/release/solve_wmec_fragments}"
if [[ -x "$RUST_SOLVER_BIN" ]]; then
  RUST_SOLVER_CMD=("$RUST_SOLVER_BIN")
else
  RUST_SOLVER_CMD=(cargo run --quiet --release --bin solve_wmec_fragments --)
fi

echo "Phase6C solver backend: rust" >> "$ROOT/solve_phase6c.log"
"${RUST_SOLVER_CMD[@]}" "${SOLVER_ARGS[@]}" >> "$ROOT/solve_phase6c.log" 2>&1

tail -n 60 "$ROOT/solve_phase6c.log"

echo
echo "===== ${LABEL}: evaluate Phase6C ====="
python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase6c.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase6c" \
  > "$ROOT/eval_phase6c.log" 2>&1

cat "$ROOT/truth_eval_phase6c.metrics.json"
record_timing "AWPhase_Phase6C_Rust_WMEC" "subset_fragments_solve_eval" "$RUN_STARTED_EPOCH"

echo
echo "===== ${LABEL}: Phase7A panel fill balanced ====="
PYTHONPATH=python python python/awphase_py/phase7a_panel_ld_fill_v2.py \
  --panel-bcf "$PANEL" \
  --backbone-local-calls-tsv "$ROOT/local_calls.phase6c.tsv" \
  --variant-json "$VARIANT_JSON" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --max-anchor-dist-bp 20000 \
  --max-anchors-each-side 4 \
  --min-panel-samples 20 \
  --min-confidence 0.85 \
  --min-margin 5 \
  --min-best-vs-second-margin 2 \
  --out-tsv "$ROOT/local_calls.phase7a.tsv" \
  --out-fill-tsv "$ROOT/fill_candidates.phase7a.tsv" \
  --out-summary-json "$ROOT/fill.phase7a.summary.json" \
  > "$ROOT/fill_phase7a.log" 2>&1

tail -n 80 "$ROOT/fill_phase7a.log"

echo
echo "===== ${LABEL}: evaluate Phase7A ====="
python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$ROOT/local_calls.phase7a.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "$ROOT/truth_eval_phase7a" \
  > "$ROOT/eval_phase7a.log" 2>&1

cat "$ROOT/truth_eval_phase7a.metrics.json"
record_timing "AWPhase_Phase7A_current" "phase6c_plus_panel_fill_eval" "$RUN_STARTED_EPOCH"

echo
echo "===== ${LABEL}: orientation-aware site comparison ====="
PYTHONPATH=python python python/awphase_py/orient_site_comparison_v1.py \
  --site-comparison-tsv "$ROOT/truth_eval_phase7a.site_comparison.tsv" \
  --out-tsv "$ROOT/truth_eval_phase7a.site_comparison.oriented.tsv" \
  --out-summary-tsv "$ROOT/truth_eval_phase7a.site_comparison.oriented.summary.tsv" \
  --min-orientation-sites 2 \
  > "$ROOT/orient_site_comparison.log" 2>&1

cat "$ROOT/truth_eval_phase7a.site_comparison.oriented.summary.tsv"

echo
echo "===== ${LABEL}: panel-fill audit ====="
PYTHONPATH=python python python/awphase_py/audit_phase7a_filled_sites_v2.py \
  --local-calls-tsv "$ROOT/local_calls.phase7a.tsv" \
  --fill-candidates-tsv "$ROOT/fill_candidates.phase7a.tsv" \
  --site-comparison-tsv "$ROOT/truth_eval_phase7a.site_comparison.tsv" \
  --out-tsv "$ROOT/panel_fill_audit.v2.tsv" \
  --out-summary-tsv "$ROOT/panel_fill_audit.v2.summary.tsv" \
  --min-orientation-sites 2 \
  > "$ROOT/audit_phase7a.log" 2>&1

cat "$ROOT/panel_fill_audit.v2.summary.tsv"

echo
echo "===== ${LABEL}: XGBoost training table ====="
PYTHONPATH=python python python/awphase_py/make_phase7a_xgboost_training_table_v1.py \
  --audit-tsv "$ROOT/panel_fill_audit.v2.tsv" \
  --out-tsv "$ROOT/phase7a_xgboost_training.tsv" \
  > "$ROOT/make_xgb_training.log" 2>&1

cat "$ROOT/make_xgb_training.log"

echo
echo "DONE: $ROOT"
