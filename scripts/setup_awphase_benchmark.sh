#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-6}"
CURL_RETRIES="${CURL_RETRIES:-5}"
PIXI_BIN="${PIXI_BIN:-$HOME/.pixi/bin/pixi}"
CHROMS=(chr1 chr20 chr22)

GIAB_TRUTH_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q"
GIAB_REF_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38"
GNOMAD_PANEL_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
GNOMAD_META_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"
BEAGLE_MAP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
ILLUMINA_BAM_URL="https://storage.googleapis.com/deepvariant/vg-case-study/HG002.novaseq.pcr-free.35x.vg-1.55.0.bam"
PACBIO_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031"
PACBIO_BAM_URL="${PACBIO_BASE}/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"

REF_GZ="references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
REF_FASTA="references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
MANIFEST="results/phase8f_manifests/phase8f_windows.window_beds.tsv"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/setup_awphase_benchmark.sh setup
  bash scripts/setup_awphase_benchmark.sh validate
  bash scripts/setup_awphase_benchmark.sh smoke
  bash scripts/setup_awphase_benchmark.sh all

Modes:
  setup     install/use Pixi, download/stage chr1/chr20/chr22 data, derive local inputs
  validate  check staged files and AWPhase benchmark bundle wiring
  smoke     run chr22:20-25 Mb AWPhase + WhatsHap Illumina + WhatsHap PacBio HiFi
  all       setup, validate, then smoke

Environment:
  THREADS=6          samtools/bcftools thread count
  PIXI_BIN=...       explicit Pixi binary path
  CURL_RETRIES=5     retry count for HTTP downloads
USAGE
}

log() {
  printf '\n===== %s =====\n' "$*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

ensure_pixi() {
  if [[ -x "$PIXI_BIN" ]]; then
    :
  elif command -v pixi >/dev/null 2>&1; then
    PIXI_BIN="$(command -v pixi)"
  else
    log "Install Pixi"
    curl -fsSL https://pixi.sh/install.sh | bash
    PIXI_BIN="$HOME/.pixi/bin/pixi"
  fi

  [[ -x "$PIXI_BIN" ]] || die "Pixi is not executable at $PIXI_BIN"

  log "Install Pixi environment"
  "$PIXI_BIN" install
}

pixi_run() {
  "$PIXI_BIN" run "$@"
}

download_url() {
  local url="$1"
  local out="$2"
  local tmp="${out}.part"

  mkdir -p "$(dirname "$out")"

  if [[ -s "$out" ]]; then
    echo "already present: $out"
    return
  fi

  echo "download: $url"
  echo "      to: $out"
  curl -fL --retry "$CURL_RETRIES" --retry-delay 5 --continue-at - --output "$tmp" "$url"
  mv "$tmp" "$out"
}

chrom_len() {
  case "$1" in
    chr1) echo 248956422 ;;
    chr20) echo 64444167 ;;
    chr22) echo 50818468 ;;
    *) die "unsupported chromosome: $1" ;;
  esac
}

illumina_out_for_chrom() {
  case "$1" in
    chr1) echo "data/raw/hg002_chr1/illumina_chrom/HG002.illumina.30x.from35x.chr1.bam" ;;
    chr20) echo "baselines/illumina_30x/HG002.illumina.chr20.30x.bam" ;;
    chr22) echo "data/raw/hg002_chr22/illumina_chrom/HG002.illumina.30x.from35x.chr22.bam" ;;
    *) die "unsupported chromosome: $1" ;;
  esac
}

pacbio_out_for_chrom() {
  local chrom="$1"
  echo "data/raw/hg002_${chrom}/pacbio_hifi/HG002.pacbio-revio-hifi.30x.${chrom}.bam"
}

truth_chr_vcf() {
  local chrom="$1"
  echo "data/truth/hg002_${chrom}/HG002_GRCh38_v5.0q_smvar.${chrom}.vcf.gz"
}

truth_chr_bed() {
  local chrom="$1"
  echo "data/truth/hg002_${chrom}/HG002_GRCh38_v5.0q_smvar.${chrom}.benchmark.bed"
}

variant_json_for_chrom() {
  local chrom="$1"
  echo "data/derived/full_chrom_variants/hg002_${chrom}_full.v5q.variants.real.json"
}

download_panels() {
  log "Download HGDP + 1KGP SHAPEIT5 panels"
  mkdir -p data/panels/hgdp_1kg_hg38/full data/panels/hgdp_1kg_hg38/meta
  for chrom in "${CHROMS[@]}"; do
    local base="hgdp1kgp_${chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
    download_url "${GNOMAD_PANEL_BASE}/${base}" "data/panels/hgdp_1kg_hg38/full/${base}"
    download_url "${GNOMAD_PANEL_BASE}/${base}.csi" "data/panels/hgdp_1kg_hg38/full/${base}.csi"
  done
  download_url "$GNOMAD_META_URL" "data/panels/hgdp_1kg_hg38/meta/gnomad_meta_updated.tsv"
}

download_truth() {
  log "Download GIAB HG002 v5.0q GRCh38 truth"
  local outdir="data/truth/hg002_v5.0q"
  mkdir -p "$outdir"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.vcf.gz" "${outdir}/HG002_GRCh38_v5.0q_smvar.vcf.gz"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi" "${outdir}/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.benchmark.bed" "${outdir}/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
  download_url "${GIAB_TRUTH_BASE}/NIST_HG002_v5.0q_variant-benchmarksets_README.md" "${outdir}/NIST_HG002_v5.0q_variant-benchmarksets_README.md"
  download_url "${GIAB_TRUTH_BASE}/checksum.md5" "${outdir}/checksum.md5"
}

stage_reference() {
  log "Download and index GRCh38 no-alt reference"
  mkdir -p references
  download_url "${GIAB_REF_BASE}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz" "$REF_GZ"

  if [[ ! -s "$REF_FASTA" ]]; then
    gzip -cd "$REF_GZ" > "${REF_FASTA}.part"
    mv "${REF_FASTA}.part" "$REF_FASTA"
  else
    echo "already present: $REF_FASTA"
  fi

  if [[ ! -s "${REF_FASTA}.fai" ]]; then
    pixi_run samtools faidx "$REF_FASTA"
  fi
}

stage_maps() {
  log "Download Beagle GRCh38 genetic maps"
  local zip="data/maps/beagle_grch38/plink.GRCh38.map.zip"
  local mapdir="data/maps/beagle_grch38/no_chr_in_chrom_field"
  mkdir -p "$mapdir"
  download_url "$BEAGLE_MAP_URL" "$zip"

  for chrom in "${CHROMS[@]}"; do
    local map="plink.${chrom}.GRCh38.map"
    if [[ -s "${mapdir}/${map}" ]]; then
      echo "already present: ${mapdir}/${map}"
    else
      unzip -j "$zip" "no_chr_in_chrom_field/${map}" -d "$mapdir"
    fi
  done
}

download_bam_slice() {
  local url="$1"
  local chrom="$2"
  local out="$3"
  local seed_fraction="$4"

  mkdir -p "$(dirname "$out")"

  if [[ -s "$out" && -s "${out}.bai" ]]; then
    echo "already present: $out"
    return
  fi

  rm -f "$out" "${out}.bai"
  pixi_run samtools view -@ "$THREADS" -b -s "$seed_fraction" "$url" "$chrom" -o "$out"
  pixi_run samtools index -@ "$THREADS" "$out"
}

stage_reads() {
  log "Download chromosome-sliced HG002 Illumina 151bp and PacBio Revio HiFi"
  for chrom in "${CHROMS[@]}"; do
    download_bam_slice "$ILLUMINA_BAM_URL" "$chrom" "$(illumina_out_for_chrom "$chrom")" "13700.857143"
    download_bam_slice "$PACBIO_BAM_URL" "$chrom" "$(pacbio_out_for_chrom "$chrom")" "13700.625"
  done

  local pacbio_meta="data/raw/hg002_chr1/pacbio_hifi"
  download_url "${PACBIO_BASE}/README_HG002-PacBio-Revio.md" "${pacbio_meta}/README_HG002-PacBio-Revio.md"
  download_url "${PACBIO_BASE}/checksums.md5" "${pacbio_meta}/checksums.md5"
}

stage_truth_chromosomes() {
  log "Split truth VCF/BED by chromosome"
  local full_vcf="data/truth/hg002_v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz"
  local full_bed="data/truth/hg002_v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed"

  for chrom in "${CHROMS[@]}"; do
    local outdir="data/truth/hg002_${chrom}"
    local vcf
    local bed
    vcf="$(truth_chr_vcf "$chrom")"
    bed="$(truth_chr_bed "$chrom")"
    mkdir -p "$outdir"

    if [[ ! -s "$vcf" ]]; then
      pixi_run bcftools view -r "$chrom" -Oz -o "$vcf" "$full_vcf"
    else
      echo "already present: $vcf"
    fi
    if [[ ! -s "${vcf}.tbi" ]]; then
      pixi_run tabix -f -p vcf "$vcf"
    fi

    if [[ ! -s "$bed" ]]; then
      pixi_run python python/awphase_py/subset_bed_window_v1.py \
        --in-bed "$full_bed" \
        --chrom "$chrom" \
        --start 1 \
        --end "$(chrom_len "$chrom")" \
        --out-bed "$bed"
    else
      echo "already present: $bed"
    fi
  done
}

derive_variant_jsons() {
  log "Build full-chromosome AWPhase variant JSONs"
  mkdir -p data/derived/full_chrom_variants
  for chrom in "${CHROMS[@]}"; do
    local out_json
    out_json="$(variant_json_for_chrom "$chrom")"
    if [[ -s "$out_json" ]]; then
      echo "already present: $out_json"
      continue
    fi

    pixi_run python python/awphase_py/vcf_window_to_variant_json_v1.py \
      --truth-vcf "$(truth_chr_vcf "$chrom")" \
      --chrom "$chrom" \
      --start 1 \
      --end "$(chrom_len "$chrom")" \
      --out-json "$out_json"
  done
}

derive_window_beds() {
  log "Build benchmark window BEDs from manifest"
  [[ -s "$MANIFEST" ]] || die "missing manifest: $MANIFEST"

  awk -F '\t' 'NR > 1 {print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $13}' "$MANIFEST" |
    while IFS=$'\t' read -r chrom start end label out_bed; do
      [[ -n "$chrom" && -n "$out_bed" ]] || continue
      if [[ -s "$out_bed" ]]; then
        echo "already present: $out_bed"
        continue
      fi
      pixi_run python python/awphase_py/subset_bed_window_v1.py \
        --in-bed "$(truth_chr_bed "$chrom")" \
        --chrom "$chrom" \
        --start "$start" \
        --end "$end" \
        --out-bed "$out_bed"
      echo "created ${label}: $out_bed"
    done
}

validate_inputs() {
  log "Validate benchmark inputs"
  ensure_pixi
  pixi_run validate-benchmark-bundles
  pixi_run phase9d-window check chr22_20_25mb

  pixi_run samtools quickcheck -v \
    "$(illumina_out_for_chrom chr1)" \
    "$(illumina_out_for_chrom chr20)" \
    "$(illumina_out_for_chrom chr22)" \
    "$(pacbio_out_for_chrom chr1)" \
    "$(pacbio_out_for_chrom chr20)" \
    "$(pacbio_out_for_chrom chr22)"

  for chrom in "${CHROMS[@]}"; do
    pixi_run bcftools index -n "data/panels/hgdp_1kg_hg38/full/hgdp1kgp_${chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf" >/dev/null
    pixi_run bcftools index -n "$(truth_chr_vcf "$chrom")" >/dev/null
  done

  echo
  du -sh data baselines references 2>/dev/null || true
  df -h .
}

setup_all_inputs() {
  ensure_pixi
  download_panels
  download_truth
  stage_reference
  stage_maps
  stage_reads
  stage_truth_chromosomes
  derive_variant_jsons
  derive_window_beds
}

run_smoke() {
  log "Run chr22_20_25mb smoke benchmark"
  ensure_pixi
  pixi_run bash scripts/run_window_awphase_and_whatshap.sh \
    chr22 20000000 24999999 chr22_20_25mb \
    "$(illumina_out_for_chrom chr22)" \
    "$(pacbio_out_for_chrom chr22)" \
    "$REF_FASTA"
}

main() {
  local mode="${1:-setup}"
  case "$mode" in
    setup)
      setup_all_inputs
      ;;
    validate)
      validate_inputs
      ;;
    smoke)
      run_smoke
      ;;
    all)
      setup_all_inputs
      validate_inputs
      run_smoke
      ;;
    -h|--help|help)
      usage
      ;;
    *)
      usage >&2
      die "unknown mode: $mode"
      ;;
  esac
}

main "$@"
