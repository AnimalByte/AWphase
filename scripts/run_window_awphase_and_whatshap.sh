#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
  echo "Usage:"
  echo "  bash scripts/run_window_awphase_and_whatshap.sh CHROM START END LABEL ILLUMINA30_BAM PACBIO_BAM_OR_NONE REF_FASTA"
  exit 1
fi

CHROM="$1"
START="$2"
END="$3"
LABEL="$4"
ILLUMINA30_BAM="$5"
PACBIO_BAM="$6"
REF_FASTA="$7"

PANEL="data/panels/hgdp_1kg_hg38/full/hgdp1kgp_${CHROM}.filtered.SNV_INDEL.phased.shapeit5.bcf"
FULL_JSON="data/derived/full_chrom_variants/hg002_${CHROM}_full.v5q.variants.real.json"
TRUTH_VCF="data/truth/hg002_${CHROM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz"
FULL_BED="data/truth/hg002_${CHROM}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed"
WIN_BED="data/truth/hg002_${CHROM}/HG002_GRCh38_v5.0q_smvar.${LABEL}.benchmark.bed"

mkdir -p results/experiments

echo "===== INPUTS ====="
echo "CHROM=$CHROM"
echo "START=$START"
echo "END=$END"
echo "LABEL=$LABEL"
echo "ILLUMINA30_BAM=$ILLUMINA30_BAM"
echo "PACBIO_BAM=$PACBIO_BAM"
echo "REF_FASTA=$REF_FASTA"
echo "PANEL=$PANEL"
echo "FULL_JSON=$FULL_JSON"
echo "TRUTH_VCF=$TRUTH_VCF"
echo "WIN_BED=$WIN_BED"

for f in "$ILLUMINA30_BAM" "$REF_FASTA" "$PANEL" "$FULL_JSON" "$TRUTH_VCF" "$FULL_BED"; do
  if [[ ! -s "$f" ]]; then
    echo "MISSING required file: $f"
    exit 1
  fi
done

echo
echo "===== subset BED ====="
PYTHONPATH=python python python/awphase_py/subset_bed_window_v1.py \
  --in-bed "$FULL_BED" \
  --chrom "$CHROM" \
  --start "$START" \
  --end "$END" \
  --out-bed "$WIN_BED"

echo
echo "===== run AWPhase Phase6C/7A ====="
AW_LABEL="${LABEL}_illumina30x_phase7a"
AW_ROOT="results/phase7a_windows/${AW_LABEL}"

bash scripts/run_phase7a_window.sh \
  "$CHROM" \
  "$START" \
  "$END" \
  "$AW_LABEL" \
  "$ILLUMINA30_BAM" \
  "$PANEL" \
  "$FULL_JSON" \
  "$TRUTH_VCF" \
  "$WIN_BED" \
  > "results/experiments/phase7a_${AW_LABEL}.log" 2>&1

tail -n 120 "results/experiments/phase7a_${AW_LABEL}.log"

echo
echo "===== make Phase7A candidate XGBoost table ====="
PYTHONPATH=python python python/awphase_py/make_phase7a_candidate_training_table_v2.py \
  --fill-candidates-tsv "$AW_ROOT/fill_candidates.phase7a.tsv" \
  --oriented-site-comparison-tsv "$AW_ROOT/truth_eval_phase7a.site_comparison.oriented.tsv" \
  --source-window "$(basename "$AW_ROOT")" \
  --out-tsv "$AW_ROOT/phase7a_candidate_xgboost_training.v2.tsv"

echo
echo "===== WhatsHap Illumina 30x baseline ====="
bash scripts/run_whatshap_baseline_window.sh \
  "$CHROM" \
  "$START" \
  "$END" \
  "${LABEL}_whatshap_illumina30x" \
  "$ILLUMINA30_BAM" \
  "$AW_ROOT/variants.window.json" \
  "$TRUTH_VCF" \
  "$WIN_BED" \
  "$REF_FASTA" \
  > "results/experiments/whatshap_${LABEL}_illumina30x.log" 2>&1

tail -n 80 "results/experiments/whatshap_${LABEL}_illumina30x.log"

if [[ "$PACBIO_BAM" != "NONE" && -s "$PACBIO_BAM" ]]; then
  echo
  echo "===== WhatsHap PacBio baseline ====="
  bash scripts/run_whatshap_baseline_window.sh \
    "$CHROM" \
    "$START" \
    "$END" \
    "${LABEL}_whatshap_pacbio" \
    "$PACBIO_BAM" \
    "$AW_ROOT/variants.window.json" \
    "$TRUTH_VCF" \
    "$WIN_BED" \
    "$REF_FASTA" \
    > "results/experiments/whatshap_${LABEL}_pacbio.log" 2>&1

  tail -n 80 "results/experiments/whatshap_${LABEL}_pacbio.log"
else
  echo
  echo "Skipping PacBio WhatsHap: PACBIO_BAM is NONE or missing."
fi

echo
echo "DONE: $LABEL"
