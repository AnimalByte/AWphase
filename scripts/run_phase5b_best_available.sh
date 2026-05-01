#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

CHROM="${1:-chr20}"
NUM="${CHROM#chr}"

OUTDIR="results/phase5b/${CHROM}_best_available_graft"
mkdir -p "$OUTDIR" results/experiments

pick_first_existing() {
  for p in "$@"; do
    if [[ -s "$p" ]]; then
      echo "$p"
      return 0
    fi
  done
  return 1
}

TAGGED="$(pick_first_existing \
  "results/experiments/${CHROM}_refaware_v3_tagged_on/local_calls.tsv" \
  "results/experiments/${CHROM}_refaware_v3_tagged/local_calls.tsv" \
  "results/experiments/${CHROM}_tagged_on/local_calls.tsv" \
  || true)"

FILL="$(pick_first_existing \
  "results/phase5/${CHROM}_tagged_spine_readbridge/sweep/d100000_p1_m0.5.local_calls.tsv" \
  "results/phase4e/sweep/${CHROM}_p4e_sing_r005_o100_b2_f1.local_calls.tsv" \
  "results/phase4e/${CHROM}_bridge_gated/local_calls.phase4e.singletons.tsv" \
  "results/phase4c/${CHROM}_wmec/local_calls.phase4c.tsv" \
  "results/phase4b/${CHROM}_wmec/local_calls.phase4b.tsv" \
  "results/phase4/${CHROM}_mec/local_calls.phase4.tsv" \
  "results/experiments/${CHROM}_refaware_v3_superreads_on/local_calls.tsv" \
  "results/experiments/${CHROM}_refaware_v3_selected_on/local_calls.tsv" \
  "results/experiments/${CHROM}_refaware_v3_on/local_calls.tsv" \
  || true)"

VARIANT_JSON="$(pick_first_existing \
  "data/derived/hg002_${NUM}_variants.real.json" \
  "data/derived/hg002_${CHROM}_variants.real.json" \
  "data/derived/HG002_${CHROM}_variants.real.json" \
  "data/derived/hg002_${NUM}_variants.json" \
  || true)"

TRUTH_VCF="$(pick_first_existing \
  "data/truth/hg002_${NUM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz" \
  "data/truth/hg002_${CHROM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz" \
  "data/truth/hg002_${NUM}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}.vcf.gz" \
  "data/truth/hg002_${CHROM}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}.vcf.gz" \
  || true)"

TRUTH_BED="$(pick_first_existing \
  "data/truth/hg002_${NUM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed" \
  "data/truth/hg002_${CHROM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed" \
  "data/truth/hg002_${NUM}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}.bed" \
  "data/truth/hg002_${CHROM}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}.bed" \
  || true)"

echo "CHROM=$CHROM"
echo "TAGGED=${TAGGED:-MISSING}"
echo "FILL=${FILL:-MISSING}"
echo "VARIANT_JSON=${VARIANT_JSON:-MISSING}"
echo "TRUTH_VCF=${TRUTH_VCF:-MISSING}"
echo "TRUTH_BED=${TRUTH_BED:-MISSING}"

missing=0
for label in TAGGED FILL VARIANT_JSON TRUTH_VCF TRUTH_BED; do
  val="${!label:-}"
  if [[ -z "$val" || ! -s "$val" ]]; then
    echo "Missing required input: $label" >&2
    missing=1
  fi
done

if [[ "$missing" = "1" ]]; then
  echo
  echo "Search hints:"
  echo "find results -path '*${CHROM}*' -name 'local_calls*.tsv' | sort"
  echo "find data/derived data/truth -path '*${NUM}*' -o -path '*${CHROM}*' | sort"
  exit 1
fi

OUT_CALLS="${OUTDIR}/local_calls.nearest_block.tsv"
GRAFT_SUMMARY="${OUTDIR}/nearest_block.summary.json"
AUDIT_TSV="${OUTDIR}/grafted_site_audit.tsv"
AUDIT_JSON="${OUTDIR}/grafted_site_audit.summary.json"
EVAL_PREFIX="${OUTDIR}/truth_eval_nearest_block_v5q"

PYTHONPATH=python python python/awphase_py/graft_filled_calls_to_tagged_blocks_v1.py \
  --tagged-local-calls-tsv "$TAGGED" \
  --filled-local-calls-tsv "$FILL" \
  --mode nearest_block \
  --max-anchor-dist-bp "${MAX_ANCHOR_DIST_BP:-100000}" \
  --out-tsv "$OUT_CALLS" \
  --out-summary-json "$GRAFT_SUMMARY"

PYTHONPATH=python python python/awphase_py/audit_phase5b_grafts_v1.py \
  --tagged-local-calls-tsv "$TAGGED" \
  --grafted-local-calls-tsv "$OUT_CALLS" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --out-tsv "$AUDIT_TSV" \
  --out-summary-json "$AUDIT_JSON"

python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$OUT_CALLS" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "$EVAL_PREFIX"

echo
echo "===== graft summary ====="
cat "$GRAFT_SUMMARY"
echo
echo "===== audit summary ====="
cat "$AUDIT_JSON"
echo
echo "===== eval metrics ====="
cat "${EVAL_PREFIX}.metrics.json"
