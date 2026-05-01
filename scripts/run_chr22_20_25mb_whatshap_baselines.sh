#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

SAMPLE="${SAMPLE:-HG002}"
REGION="${REGION:-chr22:20000000-25000000}"

VCF_IN="${VCF_IN:-data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.vcf.gz}"
VARIANT_JSON="${VARIANT_JSON:-data/derived/hg002_chr22_variants.real.json}"
TRUTH_VCF="${TRUTH_VCF:-data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.vcf.gz}"
TRUTH_BED="${TRUTH_BED:-data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.benchmark.bed}"
REF="${REF:-references/grch38_noalt_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"

OUTROOT="baselines/chr22_20_25mb"
mkdir -p "$OUTROOT"/illumina_30x "$OUTROOT"/pacbio_30x results/comparisons/chr22_20_25mb_awphase_vs_real_whatshap

echo "SAMPLE=$SAMPLE"
echo "REGION=$REGION"
echo "VCF_IN=$VCF_IN"
echo "VARIANT_JSON=$VARIANT_JSON"
echo "TRUTH_VCF=$TRUTH_VCF"
echo "TRUTH_BED=$TRUTH_BED"
echo "REF=$REF"

for p in "$VCF_IN" "$VARIANT_JSON" "$TRUTH_VCF" "$TRUTH_BED" "$REF"; do
  if [[ ! -s "$p" ]]; then
    echo "Missing required file: $p" >&2
    exit 1
  fi
done

if [[ ! -s "${REF}.fai" ]]; then
  echo "Indexing FASTA..."
  samtools faidx "$REF"
fi

if ! command -v whatshap >/dev/null 2>&1; then
  echo "Missing whatshap in current environment." >&2
  echo "Install with one of:" >&2
  echo "  conda install -c bioconda whatshap" >&2
  echo "  mamba install -c bioconda whatshap" >&2
  exit 2
fi

if ! command -v bcftools >/dev/null 2>&1; then
  echo "Missing bcftools in current environment." >&2
  exit 2
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "Missing samtools in current environment." >&2
  exit 2
fi

pick_first_existing() {
  for p in "$@"; do
    if [[ -s "$p" ]]; then
      echo "$p"
      return 0
    fi
  done
  return 1
}

discover_read_file() {
  label="$1"

  case "$label" in
    illumina)
      find data raw reads baselines -type f \
        \( -name "*.bam" -o -name "*.cram" \) 2>/dev/null \
        | grep -Ei 'chr22|HG002|hg002' \
        | grep -Ei 'illumina|novaseq|short|sr|30x|pcr-free' \
        | sort \
        | head -n 1
      ;;
    pacbio)
      find data raw reads baselines -type f \
        \( -name "*.bam" -o -name "*.cram" \) 2>/dev/null \
        | grep -Ei 'chr22|HG002|hg002' \
        | grep -Ei 'pacbio|hifi|ccs|pb' \
        | sort \
        | head -n 1
      ;;
  esac
}

ILLUMINA_BAM="${ILLUMINA_BAM:-$(discover_read_file illumina || true)}"
PACBIO_BAM="${PACBIO_BAM:-$(discover_read_file pacbio || true)}"

echo "ILLUMINA_BAM=${ILLUMINA_BAM:-MISSING}"
echo "PACBIO_BAM=${PACBIO_BAM:-MISSING}"

validate_reads_for_chr22() {
  label="$1"
  bam="$2"

  if [[ -z "${bam:-}" || ! -s "$bam" ]]; then
    echo "$label reads missing."
    return 1
  fi

  # Hard stop if the filename is clearly chr20. This is what broke the previous run.
  if echo "$bam" | grep -Eq 'chr20|chrom20'; then
    echo "ERROR: $label reads path looks like chr20, not chr22: $bam" >&2
    return 1
  fi

  # Check whether the alignment file actually has chr22 records.
  if samtools idxstats "$bam" 2>/dev/null | awk '$1=="chr22" && $3>0 {found=1} END{exit(found?0:1)}'; then
    echo "$label has chr22 mapped reads: $bam"
    return 0
  fi

  echo "ERROR: $label has no mapped chr22 reads according to samtools idxstats: $bam" >&2
  return 1
}


if ! validate_reads_for_chr22 "Illumina" "${ILLUMINA_BAM:-}"; then
  echo "Could not validate Illumina chr22 BAM/CRAM." >&2
  echo "Candidates:" >&2
  find data raw reads baselines -type f \( -name "*.bam" -o -name "*.cram" \) 2>/dev/null | grep -Ei 'HG002|hg002|chr22|illumina|novaseq|30x|pcr-free' | sort | sed -n '1,150p' >&2
  ILLUMINA_BAM=""
fi

if ! validate_reads_for_chr22 "PacBio/HiFi" "${PACBIO_BAM:-}"; then
  echo "Could not validate PacBio/HiFi chr22 BAM/CRAM." >&2
  echo "Candidates:" >&2
  find data raw reads baselines -type f \( -name "*.bam" -o -name "*.cram" \) 2>/dev/null | grep -Ei 'HG002|hg002|chr22|hifi|pacbio|ccs|pb' | sort | sed -n '1,150p' >&2
  PACBIO_BAM=""
fi

# Build unphased VCF for WhatsHap input.
UNPHASED_VCF="$OUTROOT/HG002.chr22_20_25mb.unphased_for_whatshap.vcf"
UNPHASED_VCFGZ="$OUTROOT/HG002.chr22_20_25mb.unphased_for_whatshap.vcf.gz"

PYTHONPATH=python python python/awphase_py/unphase_vcf_for_whatshap_v1.py \
  --in-vcf "$VCF_IN" \
  --out-vcf "$UNPHASED_VCF" \
  --sample "$SAMPLE"

bcftools view -Oz -o "$UNPHASED_VCFGZ" "$UNPHASED_VCF"
bcftools index -t "$UNPHASED_VCFGZ"

index_reads_if_needed() {
  bam="$1"
  if [[ "$bam" == *.bam ]]; then
    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
      echo "Indexing BAM: $bam"
      samtools index "$bam"
    fi
  elif [[ "$bam" == *.cram ]]; then
    if [[ ! -s "${bam}.crai" ]]; then
      echo "Indexing CRAM: $bam"
      samtools index "$bam"
    fi
  fi
}

run_whatshap_baseline() {
  label="$1"
  reads="$2"
  outdir="$OUTROOT/$label"

  if [[ -z "${reads:-}" || ! -s "$reads" ]]; then
    echo "Skipping $label because reads file is missing."
    return 0
  fi

  mkdir -p "$outdir"

  index_reads_if_needed "$reads"

  raw_vcf="$outdir/HG002.chr22_20_25mb.whatshap_${label}.vcf"
  phased_vcfgz="$outdir/HG002.chr22_20_25mb.whatshap_${label}.vcf.gz"
  local_calls="$outdir/local_calls.tsv"
  eval_prefix="$outdir/truth_eval"

  echo
  echo "===== Running WhatsHap $label ====="
  echo "reads=$reads"

  /usr/bin/time -v whatshap phase \
    --reference "$REF" \
    --sample "$SAMPLE" \
    --ignore-read-groups \
    --output "$raw_vcf" \
    "$UNPHASED_VCFGZ" \
    "$reads" \
    > "$outdir/whatshap.stdout.log" \
    2> "$outdir/whatshap.stderr.log"

  rc=$?
  echo "whatshap rc=$rc"
  tail -n 80 "$outdir/whatshap.stderr.log" || true

  if [[ "$rc" != "0" || ! -s "$raw_vcf" ]]; then
    echo "WhatsHap failed for $label"
    return 0
  fi

  bcftools view -Oz -o "$phased_vcfgz" "$raw_vcf"
  bcftools index -t "$phased_vcfgz"

  PYTHONPATH=python python python/awphase_py/whatshap_vcf_to_local_calls_v1.py \
    --phased-vcf "$phased_vcfgz" \
    --variant-json "$VARIANT_JSON" \
    --out-tsv "$local_calls" \
    --sample "$SAMPLE" \
    > "$outdir/convert_to_local_calls.log" 2>&1

  cat "$outdir/convert_to_local_calls.log"

  python python/awphase_py/evaluate_phase_truth.py \
    --pred-tsv "$local_calls" \
    --variant-json "$VARIANT_JSON" \
    --truth-vcf "$TRUTH_VCF" \
    --truth-bed "$TRUTH_BED" \
    --phase-column local_phase_state \
    --out-prefix "$eval_prefix" \
    > "$outdir/evaluate.log" 2>&1

  cat "${eval_prefix}.metrics.json"
}

run_whatshap_baseline "illumina_30x" "$ILLUMINA_BAM"
run_whatshap_baseline "pacbio_30x" "$PACBIO_BAM"

echo
echo "Done. Baselines written under: $OUTROOT"
