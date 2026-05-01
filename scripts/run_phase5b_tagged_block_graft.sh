#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

CHROM="${1:-chr20}"

mkdir -p "results/phase5b/${CHROM}_tagged_block_graft" results/experiments

case "$CHROM" in
  chr20)
    TAGGED="${TAGGED:-results/experiments/chr20_refaware_v3_tagged_on/local_calls.tsv}"
    FILL="${FILL:-results/phase5/tagged_spine_readbridge/sweep/d100000_p1_m0.5.local_calls.tsv}"
    VARIANT_JSON="${VARIANT_JSON:-data/derived/hg002_chr20_variants.real.json}"
    TRUTH_VCF="${TRUTH_VCF:-data/truth/hg002_chr20/HG002_GRCh38_v5.0q_smvar.chr20.vcf.gz}"
    TRUTH_BED="${TRUTH_BED:-data/truth/hg002_chr20/HG002_GRCh38_v5.0q_smvar.chr20.benchmark.bed}"
    ;;
  chr1)
    TAGGED="${TAGGED:-results/experiments/chr1_refaware_v3_tagged_on/local_calls.tsv}"
    FILL="${FILL:-results/phase5/chr1_tagged_spine_readbridge/sweep/d100000_p1_m0.5.local_calls.tsv}"
    VARIANT_JSON="${VARIANT_JSON:-data/derived/hg002_chr1_variants.real.json}"
    TRUTH_VCF="${TRUTH_VCF:-data/truth/hg002_chr1/HG002_GRCh38_v5.0q_smvar.chr1.vcf.gz}"
    TRUTH_BED="${TRUTH_BED:-data/truth/hg002_chr1/HG002_GRCh38_v5.0q_smvar.chr1.benchmark.bed}"
    ;;
  chr22)
    TAGGED="${TAGGED:-results/experiments/chr22_refaware_v3_tagged_on/local_calls.tsv}"
    FILL="${FILL:-results/phase5/chr22_tagged_spine_readbridge/sweep/d100000_p1_m0.5.local_calls.tsv}"
    VARIANT_JSON="${VARIANT_JSON:-data/derived/hg002_chr22_variants.real.json}"
    TRUTH_VCF="${TRUTH_VCF:-data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz}"
    TRUTH_BED="${TRUTH_BED:-data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.benchmark.bed}"
    ;;
  *)
    echo "Unsupported chromosome: $CHROM" >&2
    exit 2
    ;;
esac

OUTDIR="results/phase5b/${CHROM}_tagged_block_graft"
OUT_CALLS="${OUTDIR}/local_calls.nearest_block.tsv"
SUMMARY="${OUTDIR}/nearest_block.summary.json"
AUDIT_TSV="${OUTDIR}/grafted_site_audit.tsv"
AUDIT_JSON="${OUTDIR}/grafted_site_audit.summary.json"
EVAL_PREFIX="${OUTDIR}/truth_eval_nearest_block_v5q"

echo "CHROM=$CHROM"
echo "TAGGED=$TAGGED"
echo "FILL=$FILL"
echo "VARIANT_JSON=$VARIANT_JSON"
echo "TRUTH_VCF=$TRUTH_VCF"
echo "TRUTH_BED=$TRUTH_BED"

for p in "$TAGGED" "$FILL" "$VARIANT_JSON" "$TRUTH_VCF" "$TRUTH_BED"; do
  if [[ ! -s "$p" ]]; then
    echo "Missing required input: $p" >&2
    exit 1
  fi
done

PYTHONPATH=python python python/awphase_py/graft_filled_calls_to_tagged_blocks_v1.py \
  --tagged-local-calls-tsv "$TAGGED" \
  --filled-local-calls-tsv "$FILL" \
  --mode nearest_block \
  --max-anchor-dist-bp "${MAX_ANCHOR_DIST_BP:-100000}" \
  --out-tsv "$OUT_CALLS" \
  --out-summary-json "$SUMMARY"

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
cat "$SUMMARY"
echo
echo "===== audit summary ====="
cat "$AUDIT_JSON"
echo
echo "===== eval metrics ====="
cat "${EVAL_PREFIX}.metrics.json"
