#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

SAMPLE="HG002"
CHROM="chr22"
START=20000000
END=25000000

REF="references/grch38_noalt_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
VARIANT_JSON="data/derived/hg002_chr22_variants.real.json"
TRUTH_VCF="data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.vcf.gz"
TRUTH_BED="data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.benchmark.bed"

UNPHASED_VCFGZ="baselines/chr22_20_25mb/HG002.chr22_20_25mb.unphased_for_whatshap.vcf.gz"

OUTROOT="results/depth_curve/chr22_20_25mb"
mkdir -p "$OUTROOT"

if [[ ! -s "$UNPHASED_VCFGZ" ]]; then
  mkdir -p baselines/chr22_20_25mb
  PYTHONPATH=python python python/awphase_py/unphase_vcf_for_whatshap_v1.py \
    --in-vcf "$TRUTH_VCF" \
    --out-vcf baselines/chr22_20_25mb/HG002.chr22_20_25mb.unphased_for_whatshap.vcf \
    --sample "$SAMPLE"

  bcftools view -Oz \
    -o "$UNPHASED_VCFGZ" \
    baselines/chr22_20_25mb/HG002.chr22_20_25mb.unphased_for_whatshap.vcf

  bcftools index -t "$UNPHASED_VCFGZ"
fi

run_one_depth() {
  label="$1"
  bam="$2"

  echo
  echo "=============================="
  echo "DEPTH $label"
  echo "BAM=$bam"
  echo "=============================="

  if [[ ! -s "$bam" ]]; then
    echo "Missing BAM: $bam" >&2
    return 1
  fi

  if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
    samtools index "$bam"
  fi

  echo "idxstats:"
  samtools idxstats "$bam" | awk '$1=="chr22" || NR<=3'

  ########################
  # Phase6C
  ########################

  p6dir="$OUTROOT/phase6c_${label}"
  mkdir -p "$p6dir"

  echo
  echo "---- Phase6C $label: build fragments ----"
  /usr/bin/time -v env PYTHONPATH=python python python/awphase_py/build_template_fragments_from_bam_v1.py \
    --bam "$bam" \
    --variant-json "$VARIANT_JSON" \
    --chrom "$CHROM" \
    --start "$START" \
    --end "$END" \
    --min-mapq 20 \
    --min-baseq 15 \
    --min-sites-per-fragment 2 \
    --out-tsv "$p6dir/fragments.tsv" \
    --out-summary-json "$p6dir/fragments.summary.json" \
    > "$p6dir/build_fragments.log" 2>&1

  tail -n 60 "$p6dir/build_fragments.log"

  echo
  echo "---- Phase6C $label: solve ----"
  /usr/bin/time -v env PYTHONPATH=python python python/awphase_py/solve_wmec_fragments_v1.py \
    --fragments-tsv "$p6dir/fragments.tsv" \
    --variant-json "$VARIANT_JSON" \
    --max-exact-sites 18 \
    --max-component-sites 1024 \
    --local-refine-iters 8 \
    --out-local-calls-tsv "$p6dir/local_calls.phase6c.tsv" \
    --out-components-tsv "$p6dir/components.tsv" \
    --out-summary-json "$p6dir/solve.summary.json" \
    > "$p6dir/solve.log" 2>&1

  tail -n 80 "$p6dir/solve.log"

  echo
  echo "---- Phase6C $label: evaluate ----"
  python python/awphase_py/evaluate_phase_truth.py \
    --pred-tsv "$p6dir/local_calls.phase6c.tsv" \
    --variant-json "$VARIANT_JSON" \
    --truth-vcf "$TRUTH_VCF" \
    --truth-bed "$TRUTH_BED" \
    --phase-column local_phase_state \
    --out-prefix "$p6dir/truth_eval_phase6c" \
    > "$p6dir/evaluate.log" 2>&1

  cat "$p6dir/truth_eval_phase6c.metrics.json"

  ########################
  # WhatsHap
  ########################

  whdir="$OUTROOT/whatshap_${label}"
  mkdir -p "$whdir"

  raw_vcf="$whdir/HG002.chr22_20_25mb.whatshap_${label}.vcf"
  phased_vcfgz="$whdir/HG002.chr22_20_25mb.whatshap_${label}.vcf.gz"

  echo
  echo "---- WhatsHap $label ----"
  /usr/bin/time -v whatshap phase \
    --reference "$REF" \
    --sample "$SAMPLE" \
    --ignore-read-groups \
    --output "$raw_vcf" \
    "$UNPHASED_VCFGZ" \
    "$bam" \
    > "$whdir/whatshap.stdout.log" \
    2> "$whdir/whatshap.stderr.log"

  tail -n 80 "$whdir/whatshap.stderr.log"

  bcftools view -Oz -o "$phased_vcfgz" "$raw_vcf"
  bcftools index -t "$phased_vcfgz"

  PYTHONPATH=python python python/awphase_py/whatshap_vcf_to_local_calls_v1.py \
    --phased-vcf "$phased_vcfgz" \
    --variant-json "$VARIANT_JSON" \
    --out-tsv "$whdir/local_calls.tsv" \
    --sample "$SAMPLE" \
    > "$whdir/convert_to_local_calls.log" 2>&1

  cat "$whdir/convert_to_local_calls.log"

  python python/awphase_py/evaluate_phase_truth.py \
    --pred-tsv "$whdir/local_calls.tsv" \
    --variant-json "$VARIANT_JSON" \
    --truth-vcf "$TRUTH_VCF" \
    --truth-bed "$TRUTH_BED" \
    --phase-column local_phase_state \
    --out-prefix "$whdir/truth_eval" \
    > "$whdir/evaluate.log" 2>&1

  cat "$whdir/truth_eval.metrics.json"
}

run_one_depth "20x_from35x" "data/raw/hg002_chr22/downsampled/HG002.novaseq.pcr-free.20x.from35x.chr22_20_25mb.bam"
run_one_depth "25x_from35x" "data/raw/hg002_chr22/downsampled/HG002.novaseq.pcr-free.25x.from35x.chr22_20_25mb.bam"
run_one_depth "30x_from35x" "data/raw/hg002_chr22/downsampled/HG002.novaseq.pcr-free.30x.from35x.chr22_20_25mb.bam"

echo
echo "Depth curve complete."
