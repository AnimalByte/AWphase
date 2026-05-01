#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 6 ]; then
  echo "usage: $0 <exp_name> <config> <variant_json> <results_subdir> <truth_vcf> <truth_bed> [benchmark extra args...]"
  exit 1
fi

EXP_NAME="$1"
CONFIG="$2"
VARIANT_JSON="$3"
RESULTS_SUBDIR="$4"
TRUTH_VCF="$5"
TRUTH_BED="$6"
shift 6

OUTDIR="results/experiments/${EXP_NAME}"
RUNDIR="results/runs/HG002/${RESULTS_SUBDIR}"

mkdir -p "$OUTDIR"

./target/debug/awphase-cli --config "$CONFIG" benchmark "$@"

python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "${RUNDIR}/local_calls.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "${OUTDIR}/truth_eval_local"

cp "${RUNDIR}/local_window_states.json" "${OUTDIR}/"
cp "${RUNDIR}/local_calls.tsv" "${OUTDIR}/"
cp "${RUNDIR}/benchmark_summary.json" "${OUTDIR}/"
cp "${RUNDIR}/run_report.json" "${OUTDIR}/"
cp "${RUNDIR}/candidate_set.json" "${OUTDIR}/" || true
cp "${RUNDIR}/boundary_decisions.json" "${OUTDIR}/" || true
cp "${RUNDIR}/derived_boundary_scores.json" "${OUTDIR}/" || true
cp "${RUNDIR}/derived_block_summaries.json" "${OUTDIR}/" || true

echo "Captured ${EXP_NAME} -> ${OUTDIR}"
