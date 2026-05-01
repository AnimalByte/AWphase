#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

VARIANT_JSON="data/derived/hg002_chr22_variants.real.json"
TRUTH_VCF="data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.vcf.gz"
TRUTH_BED="data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.benchmark.bed"

OUTROOT="results/comparisons/chr22_20_25mb_awphase_vs_whatshap"
mkdir -p "$OUTROOT"

BRANCHES="$OUTROOT/branches.tsv"
SUMMARY="$OUTROOT/summary.tsv"

echo -e "name\tpath" > "$BRANCHES"

add_branch() {
  name="$1"
  path="$2"
  if [[ -s "$path" ]]; then
    echo -e "${name}\t${path}" >> "$BRANCHES"
  else
    echo "missing, skipped: $name -> $path" >&2
  fi
}

# AWPhase baselines / current candidates
add_branch "awphase_tagged" \
  "results/experiments/chr22_refaware_v3_tagged_on/local_calls.tsv"

add_branch "awphase_selected" \
  "results/experiments/chr22_refaware_v3_selected_on/local_calls.tsv"

add_branch "awphase_superreads" \
  "results/experiments/chr22_refaware_v3_superreads_on/local_calls.tsv"

add_branch "awphase_readbridge_superreads_graft_d1000" \
  "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d1000/grafted.tsv"

add_branch "awphase_readbridge_superreads_graft_d2500" \
  "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d2500/grafted.tsv"

add_branch "awphase_readbridge_superreads_graft_d5000" \
  "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d5000/grafted.tsv"

# Common possible WhatsHap names. These are guesses; missing paths are skipped.
add_branch "whatshap_illumina" \
  "results/experiments/chr22_whatshap_illumina/local_calls.tsv"

add_branch "whatshap_illumina_on" \
  "results/experiments/chr22_whatshap_illumina_on/local_calls.tsv"

add_branch "whatshap_hifi" \
  "results/experiments/chr22_whatshap_hifi/local_calls.tsv"

add_branch "whatshap_hifi_on" \
  "results/experiments/chr22_whatshap_hifi_on/local_calls.tsv"

add_branch "whatshap_pacbio" \
  "results/experiments/chr22_whatshap_pacbio/local_calls.tsv"

add_branch "whatshap_pacbio_on" \
  "results/experiments/chr22_whatshap_pacbio_on/local_calls.tsv"

# Auto-add any other local_calls.tsv file with whatshap/hifi/pacbio/illumina in path.
find results -type f -name "local_calls*.tsv" 2>/dev/null \
  | grep -Ei 'chr22.*(whatshap|whats|hifi|pacbio|illumina)' \
  | sort \
  | while read -r p; do
      safe="$(echo "$p" | sed 's#^results/##; s#/#__#g; s#[^A-Za-z0-9_.-]#_#g; s#\.tsv$##')"
      if ! grep -Fq "$p" "$BRANCHES"; then
        echo -e "${safe}\t${p}" >> "$BRANCHES"
      fi
    done

echo -e "name\tpath\tn_pred_sites_nonzero\tn_exact_overlap_sites_phased\thamming_denominator\thamming_error_rate\tswitch_denominator\tswitch_error_rate\tphased_site_accuracy_pct\ttruth_correct_pct\traw_block_n50_bp\tmax_block_span_bp\tmedian_block_span_bp" > "$SUMMARY"

tail -n +2 "$BRANCHES" | while IFS=$'\t' read -r name path; do
  safe="${name//[^A-Za-z0-9_.-]/_}"
  prefix="$OUTROOT/${safe}.eval"

  echo "Evaluating $name -> $path"

  python python/awphase_py/evaluate_phase_truth.py \
    --pred-tsv "$path" \
    --variant-json "$VARIANT_JSON" \
    --truth-vcf "$TRUTH_VCF" \
    --truth-bed "$TRUTH_BED" \
    --phase-column local_phase_state \
    --out-prefix "$prefix" \
    > "$OUTROOT/${safe}.eval.log" 2>&1

  python - <<PY >> "$SUMMARY"
import json
m=json.load(open("${prefix}.metrics.json"))
print("\\t".join(map(str, [
  "$name",
  "$path",
  m.get("n_pred_sites_nonzero"),
  m.get("n_exact_overlap_sites_phased"),
  m.get("hamming_denominator"),
  m.get("hamming_error_rate"),
  m.get("switch_denominator"),
  m.get("switch_error_rate"),
  m.get("phased_site_accuracy_pct"),
  m.get("truth_correct_pct"),
  m.get("raw_block_n50_bp"),
  m.get("max_block_span_bp"),
  m.get("median_block_span_bp"),
])))
PY
done

echo
echo "===== all branches ====="
column -t "$SUMMARY"

echo
echo "===== best hamming ====="
column -t "$SUMMARY" | sort -k6,6n | sed -n '1,30p'

echo
echo "===== best switch ====="
column -t "$SUMMARY" | sort -k8,8n | sed -n '1,30p'

echo
echo "Wrote:"
echo "  $BRANCHES"
echo "  $SUMMARY"
