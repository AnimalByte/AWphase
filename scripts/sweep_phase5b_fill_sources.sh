#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

CHROM="${1:-chr1}"
NUM="${CHROM#chr}"

OUTROOT="results/phase5b/${CHROM}_fill_source_sweep"
mkdir -p "$OUTROOT" results/experiments

pick_first_existing() {
  for p in "$@"; do
    if [[ -s "$p" ]]; then
      echo "$p"
      return 0
    fi
  done
  return 1
}

safe_name() {
  echo "$1" | sed 's#^results/##; s#/#__#g; s#[^A-Za-z0-9_.-]#_#g; s#\.tsv$##'
}

TAGGED="$(pick_first_existing \
  "results/experiments/${CHROM}_refaware_v3_tagged_on/local_calls.tsv" \
  "results/experiments/${CHROM}_refaware_v3_tagged/local_calls.tsv" \
  || true)"

VARIANT_JSON="$(pick_first_existing \
  "data/derived/hg002_${CHROM}_variants.real.json" \
  "data/derived/hg002_${NUM}_variants.real.json" \
  "data/derived/hg002_${CHROM}_variants.json" \
  "data/derived/hg002_${NUM}_variants.json" \
  || true)"

TRUTH_DIR_A="data/truth/hg002_${CHROM}"
TRUTH_DIR_B="data/truth/hg002_${NUM}"

# Strict truth VCF discovery: only search the chromosome-specific truth directory.
TRUTH_VCF="$(pick_first_existing \
  "${TRUTH_DIR_A}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz" \
  "${TRUTH_DIR_A}/HG002_GRCh38_v5.0q_smvar.${CHROM}_20_25mb.vcf.gz" \
  "${TRUTH_DIR_A}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}.vcf.gz" \
  "${TRUTH_DIR_A}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}_20_25mb.vcf.gz" \
  "${TRUTH_DIR_B}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz" \
  "${TRUTH_DIR_B}/HG002_GRCh38_v5.0q_smvar.${CHROM}_20_25mb.vcf.gz" \
  "${TRUTH_DIR_B}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}.vcf.gz" \
  "${TRUTH_DIR_B}/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.${CHROM}_20_25mb.vcf.gz" \
  || true)"

if [[ -z "${TRUTH_VCF:-}" ]]; then
  TRUTH_VCF="$(find "$TRUTH_DIR_A" "$TRUTH_DIR_B" -type f -name "*.${CHROM}*.vcf.gz" 2>/dev/null | sort | head -n 1 || true)"
fi

# Strict BED discovery: do not match chr1 files just because they contain GRCh38_1_22.
TRUTH_BED="$(pick_first_existing \
  "${TRUTH_DIR_A}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed" \
  "${TRUTH_DIR_A}/HG002_GRCh38_v5.0q_smvar.${CHROM}_20_25mb.benchmark.bed" \
  "${TRUTH_DIR_A}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}.bed" \
  "${TRUTH_DIR_A}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}_20_25mb.bed" \
  "${TRUTH_DIR_B}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed" \
  "${TRUTH_DIR_B}/HG002_GRCh38_v5.0q_smvar.${CHROM}_20_25mb.benchmark.bed" \
  "${TRUTH_DIR_B}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}.bed" \
  "${TRUTH_DIR_B}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.${CHROM}_20_25mb.bed" \
  || true)"

if [[ -z "${TRUTH_BED:-}" ]]; then
  TRUTH_BED="$(find "$TRUTH_DIR_A" "$TRUTH_DIR_B" -type f \( -name "*.${CHROM}.bed" -o -name "*.${CHROM}_*.bed" \) 2>/dev/null | grep -E 'benchmark|noinconsistent' | sort | head -n 1 || true)"
fi

OBS="$(pick_first_existing \
  "results/experiments/${CHROM}_read_frontend/readobs_refaware_v3.obs_debug.tsv" \
  "results/experiments/${CHROM}_refaware_v3_on/readobs_refaware_v3.obs_debug.tsv" \
  "results/experiments/${CHROM}_read_frontend/readobs_refaware_v3.obs_debug.tsv" \
  || true)"

echo "CHROM=$CHROM"
echo "TAGGED=${TAGGED:-MISSING}"
echo "VARIANT_JSON=${VARIANT_JSON:-MISSING}"
echo "TRUTH_VCF=${TRUTH_VCF:-MISSING}"
echo "TRUTH_BED=${TRUTH_BED:-MISSING}"
echo "OBS=${OBS:-MISSING_OPTIONAL}"

missing=0
for label in TAGGED VARIANT_JSON TRUTH_VCF TRUTH_BED; do
  val="${!label:-}"
  if [[ -z "$val" || ! -s "$val" ]]; then
    echo "Missing required input: $label=${val:-MISSING}" >&2
    missing=1
  fi
done

if [[ "$missing" = "1" ]]; then
  echo
  echo "Truth candidates for $CHROM:"
  find data/truth -type f -path "*hg002_${CHROM}*" -o -path "*hg002_${NUM}*" 2>/dev/null | sort || true
  echo
  echo "Local call candidates for $CHROM:"
  find results -path "*${CHROM}*" -name "local_calls*.tsv" | sort || true
  exit 1
fi

SUMMARY="$OUTROOT/summary.tsv"
echo -e "fill_name\tfill_path\tmode\tn_pred_sites_nonzero\tn_exact_overlap_sites_phased\thamming_denominator\thamming_error_rate\tswitch_denominator\tswitch_error_rate\tphased_site_accuracy_pct\ttruth_correct_pct\traw_block_n50_bp\tgrafted_sites\ttruth_comparable_grafted_sites\tsame_orientation_accuracy\topposite_orientation_accuracy" > "$SUMMARY"

# Baseline.
python python/awphase_py/evaluate_phase_truth.py \
  --pred-tsv "$TAGGED" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --truth-bed "$TRUTH_BED" \
  --phase-column local_phase_state \
  --out-prefix "$OUTROOT/tagged_baseline.eval" \
  > "$OUTROOT/tagged_baseline.eval.log" 2>&1

python - <<PY >> "$SUMMARY"
import json
m=json.load(open("$OUTROOT/tagged_baseline.eval.metrics.json"))
print("\\t".join(map(str, [
  "TAGGED_BASELINE",
  "$TAGGED",
  "baseline",
  m.get("n_pred_sites_nonzero"),
  m.get("n_exact_overlap_sites_phased"),
  m.get("hamming_denominator"),
  m.get("hamming_error_rate"),
  m.get("switch_denominator"),
  m.get("switch_error_rate"),
  m.get("phased_site_accuracy_pct"),
  m.get("truth_correct_pct"),
  m.get("raw_block_n50_bp"),
  "",
  "",
  "",
  "",
])))
PY

mapfile -t FILLS < <(
  find results -path "*${CHROM}*" -name "local_calls*.tsv" | sort | \
    grep -v "phase5b" | \
    grep -v "tagged_block_graft" | \
    grep -v "best_available_graft" | \
    grep -v "${CHROM}_refaware_v3_tagged_on/local_calls.tsv"
)

if [[ "${#FILLS[@]}" -eq 0 ]]; then
  echo "No fill candidates found for $CHROM" >&2
  exit 2
fi

for FILL in "${FILLS[@]}"; do
  name="$(safe_name "$FILL")"
  work="$OUTROOT/$name"
  mkdir -p "$work"

  echo
  echo "===== FILL $name ====="
  echo "FILL=$FILL"

  python python/awphase_py/evaluate_phase_truth.py \
    --pred-tsv "$FILL" \
    --variant-json "$VARIANT_JSON" \
    --truth-vcf "$TRUTH_VCF" \
    --truth-bed "$TRUTH_BED" \
    --phase-column local_phase_state \
    --out-prefix "$work/fill_alone.eval" \
    > "$work/fill_alone.eval.log" 2>&1 || {
      echo "fill eval failed for $FILL; skipping" >&2
      continue
    }

  if [[ -n "${OBS:-}" && -s "$OBS" ]]; then
    PYTHONPATH=python python python/awphase_py/build_tagged_spine_readbridge_fill_v1.py \
      --tagged-local-calls-tsv "$TAGGED" \
      --fill-local-calls-tsv "$FILL" \
      --obs-tsv "$OBS" \
      --max-anchor-dist-bp 100000 \
      --min-bridge-pairs 1 \
      --min-orientation-margin 0.5 \
      --allow-one-sided-anchors \
      --out-tsv "$work/readbridge_fill.tsv" \
      --out-summary-json "$work/readbridge_fill.summary.json" \
      > "$work/readbridge_fill.log" 2>&1 || {
        echo "readbridge fill failed for $FILL; using direct fill" >&2
        cp "$FILL" "$work/readbridge_fill.tsv"
      }
    GRAFT_INPUT="$work/readbridge_fill.tsv"
  else
    echo "No OBS file found; using direct fill without readbridge orientation"
    GRAFT_INPUT="$FILL"
  fi

  for mode in nearest_block same_block_only; do
    PYTHONPATH=python python python/awphase_py/graft_filled_calls_to_tagged_blocks_v1.py \
      --tagged-local-calls-tsv "$TAGGED" \
      --filled-local-calls-tsv "$GRAFT_INPUT" \
      --mode "$mode" \
      --max-anchor-dist-bp 100000 \
      --out-tsv "$work/grafted_${mode}.tsv" \
      --out-summary-json "$work/grafted_${mode}.summary.json" \
      > "$work/grafted_${mode}.log" 2>&1 || {
        echo "graft failed for $FILL mode=$mode; skipping" >&2
        continue
      }

    PYTHONPATH=python python python/awphase_py/audit_phase5b_grafts_v1.py \
      --tagged-local-calls-tsv "$TAGGED" \
      --grafted-local-calls-tsv "$work/grafted_${mode}.tsv" \
      --variant-json "$VARIANT_JSON" \
      --truth-vcf "$TRUTH_VCF" \
      --out-tsv "$work/grafted_${mode}.audit.tsv" \
      --out-summary-json "$work/grafted_${mode}.audit.summary.json" \
      > "$work/grafted_${mode}.audit.log" 2>&1 || true

    python python/awphase_py/evaluate_phase_truth.py \
      --pred-tsv "$work/grafted_${mode}.tsv" \
      --variant-json "$VARIANT_JSON" \
      --truth-vcf "$TRUTH_VCF" \
      --truth-bed "$TRUTH_BED" \
      --phase-column local_phase_state \
      --out-prefix "$work/grafted_${mode}.eval" \
      > "$work/grafted_${mode}.eval.log" 2>&1 || {
        echo "eval failed for $FILL mode=$mode; skipping" >&2
        continue
      }

    python - <<PY >> "$SUMMARY"
import json, pathlib
m=json.load(open("$work/grafted_${mode}.eval.metrics.json"))
g=json.load(open("$work/grafted_${mode}.summary.json"))
ap=pathlib.Path("$work/grafted_${mode}.audit.summary.json")
a=json.load(open(ap)) if ap.exists() else {}
print("\\t".join(map(str, [
  "$name",
  "$FILL",
  "$mode",
  m.get("n_pred_sites_nonzero"),
  m.get("n_exact_overlap_sites_phased"),
  m.get("hamming_denominator"),
  m.get("hamming_error_rate"),
  m.get("switch_denominator"),
  m.get("switch_error_rate"),
  m.get("phased_site_accuracy_pct"),
  m.get("truth_correct_pct"),
  m.get("raw_block_n50_bp"),
  g.get("grafted_sites"),
  a.get("truth_comparable_grafted_sites"),
  a.get("same_orientation_accuracy"),
  a.get("opposite_orientation_accuracy"),
])))
PY
  done
done

echo
echo "===== best hamming ====="
column -t "$SUMMARY" | sort -k7,7n | sed -n '1,30p'

echo
echo "===== best switch ====="
column -t "$SUMMARY" | sort -k9,9n | sed -n '1,30p'

echo
echo "Wrote: $SUMMARY"
