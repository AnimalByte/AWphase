#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-12}"
CURL_RETRIES="${CURL_RETRIES:-5}"
PIXI_BIN="${PIXI_BIN:-$HOME/.pixi/bin/pixi}"

CHROM="chr6"
CHROM_LEN="170805979"
REGION_START="25000000"
REGION_END="34999999"
REGION="${CHROM}:${REGION_START}-${REGION_END}"
SPLIT="hard_chr6"

GIAB_TRUTH_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q"
GNOMAD_PANEL_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
BEAGLE_MAP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
ILLUMINA_BAM_URL="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
PACBIO_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031"
PACBIO_BAM_URL="${PACBIO_BASE}/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"

PANEL="data/panels/hgdp_1kg_hg38/full/hgdp1kgp_${CHROM}.filtered.SNV_INDEL.phased.shapeit5.bcf"
MAP_ZIP="data/maps/beagle_grch38/plink.GRCh38.map.zip"
MAP="data/maps/beagle_grch38/no_chr_in_chrom_field/plink.${CHROM}.GRCh38.map"

FULL_TRUTH_VCF="data/truth/hg002_v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz"
FULL_TRUTH_BED="data/truth/hg002_v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
TRUTH_DIR="data/truth/hg002_${CHROM}"
TRUTH_VCF="${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.${CHROM}.vcf.gz"
TRUTH_BED="${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.${CHROM}.benchmark.bed"
VARIANT_JSON="data/derived/full_chrom_variants/hg002_${CHROM}_full.v5q.variants.real.json"

ILLUMINA_DIR="data/raw/hg002_${CHROM}/illumina_mhc"
ILLUMINA_35="${ILLUMINA_DIR}/HG002.illumina.35x.${CHROM}_25_35mb.bam"
ILLUMINA_30="${ILLUMINA_DIR}/HG002.illumina.30x.from35x.${CHROM}_25_35mb.bam"

PACBIO_DIR="data/raw/hg002_${CHROM}/pacbio_hifi"
PACBIO_48="${PACBIO_DIR}/HG002.pacbio-revio-hifi.48x.${CHROM}_25_35mb.bam"
PACBIO_30="${PACBIO_DIR}/HG002.pacbio-revio-hifi.30x.${CHROM}_25_35mb.bam"

MANIFEST="results/phase8f_manifests/chr6_mhc_windows.window_beds.tsv"

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
    die "Pixi not found. Install pixi first or set PIXI_BIN."
  fi

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

stage_panel() {
  log "Stage chr6 HGDP + 1KGP SHAPEIT5 panel"
  local base
  base="$(basename "$PANEL")"
  download_url "${GNOMAD_PANEL_BASE}/${base}" "$PANEL"
  download_url "${GNOMAD_PANEL_BASE}/${base}.csi" "${PANEL}.csi"
  pixi_run bcftools index -n "$PANEL" >/dev/null
}

stage_map() {
  log "Stage chr6 Beagle GRCh38 map"
  download_url "$BEAGLE_MAP_URL" "$MAP_ZIP"
  if [[ -s "$MAP" ]]; then
    echo "already present: $MAP"
    return
  fi

  mkdir -p "$(dirname "$MAP")"
  pixi_run python - "$MAP_ZIP" "$MAP" <<'PY'
import sys
import zipfile

zip_path, out_path = sys.argv[1], sys.argv[2]
member = "no_chr_in_chrom_field/plink.chr6.GRCh38.map"
with zipfile.ZipFile(zip_path) as zf:
    data = zf.read(member)
with open(out_path, "wb") as out:
    out.write(data)
print({"out": out_path, "bytes": len(data)})
PY
}

stage_truth() {
  log "Stage chr6 GIAB HG002 v5.0q truth"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.vcf.gz" "$FULL_TRUTH_VCF"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi" "${FULL_TRUTH_VCF}.tbi"
  download_url "${GIAB_TRUTH_BASE}/HG002_GRCh38_v5.0q_smvar.benchmark.bed" "$FULL_TRUTH_BED"

  mkdir -p "$TRUTH_DIR"
  if [[ ! -s "$TRUTH_VCF" ]]; then
    pixi_run bcftools view -r "$CHROM" -Oz -o "$TRUTH_VCF" "$FULL_TRUTH_VCF"
  else
    echo "already present: $TRUTH_VCF"
  fi
  if [[ ! -s "${TRUTH_VCF}.tbi" ]]; then
    pixi_run tabix -f -p vcf "$TRUTH_VCF"
  fi

  if [[ ! -s "$TRUTH_BED" ]]; then
    PYTHONPATH=python pixi_run python python/awphase_py/subset_bed_window_v1.py \
      --in-bed "$FULL_TRUTH_BED" \
      --chrom "$CHROM" \
      --start 1 \
      --end "$CHROM_LEN" \
      --out-bed "$TRUTH_BED"
  else
    echo "already present: $TRUTH_BED"
  fi
}

stage_variant_json() {
  log "Build chr6 AWPhase variant JSON"
  if [[ -s "$VARIANT_JSON" ]]; then
    echo "already present: $VARIANT_JSON"
    return
  fi

  mkdir -p "$(dirname "$VARIANT_JSON")"
  PYTHONPATH=python pixi_run python python/awphase_py/vcf_window_to_variant_json_v1.py \
    --truth-vcf "$TRUTH_VCF" \
    --chrom "$CHROM" \
    --start 1 \
    --end "$CHROM_LEN" \
    --out-json "$VARIANT_JSON"
}

download_slice_then_downsample() {
  local url="$1"
  local raw_bam="$2"
  local downsampled_bam="$3"
  local seed_fraction="$4"
  local label="$5"

  mkdir -p "$(dirname "$raw_bam")"

  if [[ ! -s "$raw_bam" ]]; then
    echo "download sliced ${label}: $REGION"
    pixi_run samtools view -@ "$THREADS" -b "$url" "$REGION" -o "$raw_bam"
  else
    echo "already present: $raw_bam"
  fi
  if [[ ! -s "${raw_bam}.bai" && ! -s "${raw_bam%.bam}.bai" ]]; then
    pixi_run samtools index -@ "$THREADS" "$raw_bam"
  fi

  if [[ ! -s "$downsampled_bam" ]]; then
    echo "downsample local ${label} slice to 30x: $downsampled_bam"
    pixi_run samtools view -@ "$THREADS" -b -s "$seed_fraction" "$raw_bam" -o "$downsampled_bam"
  else
    echo "already present: $downsampled_bam"
  fi
  if [[ ! -s "${downsampled_bam}.bai" && ! -s "${downsampled_bam%.bam}.bai" ]]; then
    pixi_run samtools index -@ "$THREADS" "$downsampled_bam"
  fi
}

stage_reads() {
  log "Stage chr6 MHC read slices, then downsample slices to 30x"
  download_slice_then_downsample "$ILLUMINA_BAM_URL" "$ILLUMINA_35" "$ILLUMINA_30" "42.857142" "Illumina 151bp 35x"
  download_slice_then_downsample "$PACBIO_BAM_URL" "$PACBIO_48" "$PACBIO_30" "42.625" "PacBio Revio HiFi 48x"
  download_url "${PACBIO_BASE}/README_HG002-PacBio-Revio.md" "${PACBIO_DIR}/README_HG002-PacBio-Revio.md"
  download_url "${PACBIO_BASE}/checksums.md5" "${PACBIO_DIR}/checksums.md5"
}

stage_window_beds() {
  log "Build chr6 hard-region benchmark BEDs"
  PYTHONPATH=python pixi_run python python/awphase_py/subset_bed_window_v1.py \
    --in-bed "$TRUTH_BED" \
    --chrom "$CHROM" \
    --start 25000000 \
    --end 29999999 \
    --out-bed "${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.chr6_25_30mb.benchmark.bed"

  PYTHONPATH=python pixi_run python python/awphase_py/subset_bed_window_v1.py \
    --in-bed "$TRUTH_BED" \
    --chrom "$CHROM" \
    --start 30000000 \
    --end 34999999 \
    --out-bed "${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.chr6_30_35mb.benchmark.bed"
}

write_manifest() {
  log "Write chr6 MHC manifest"
  mkdir -p "$(dirname "$MANIFEST")"
  PYTHONPATH=python pixi_run python - \
    "$MANIFEST" \
    "$ILLUMINA_30" \
    "$PANEL" \
    "$VARIANT_JSON" \
    "$TRUTH_VCF" \
    "$TRUTH_DIR" \
    "$SPLIT" <<'PY'
import csv
import sys
from pathlib import Path

import pysam

manifest, bam_path, panel, variant_json, truth_vcf, truth_dir, split = sys.argv[1:]
chrom = "chr6"
windows = [
    (25000000, 29999999, "chr6_25_30mb"),
    (30000000, 34999999, "chr6_30_35mb"),
]

def is_het_01(gt):
    return gt is not None and len(gt) == 2 and None not in gt and set(gt) == {0, 1}

def bed_intervals(path):
    n = 0
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n

with pysam.AlignmentFile(bam_path) as bam, pysam.VariantFile(truth_vcf) as vcf:
    rows = []
    for start, end, label in windows:
        truth_bed = f"{truth_dir}/HG002_GRCh38_v5.0q_smvar.{label}.benchmark.bed"
        reads = bam.count(chrom, start - 1, end)
        snps = 0
        for rec in vcf.fetch(chrom, start - 1, end):
            if rec.pos < start or rec.pos > end:
                continue
            if rec.alts is None or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue
            gt = rec.samples["HG002"].get("GT")
            if is_het_01(gt):
                snps += 1
        rows.append({
            "split": split,
            "chrom": chrom,
            "start": start,
            "end": end,
            "label": label,
            "reads": reads,
            "snp_variant_count": snps,
            "benchmark_bed_intervals": bed_intervals(truth_bed),
            "bam": bam_path,
            "panel": panel,
            "json": variant_json,
            "truth_vcf": truth_vcf,
            "truth_bed": truth_bed,
        })

fields = [
    "split",
    "chrom",
    "start",
    "end",
    "label",
    "reads",
    "snp_variant_count",
    "benchmark_bed_intervals",
    "bam",
    "panel",
    "json",
    "truth_vcf",
    "truth_bed",
]
Path(manifest).parent.mkdir(parents=True, exist_ok=True)
with open(manifest, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
print({"manifest": manifest, "rows": len(rows)})
PY
}

validate_outputs() {
  log "Validate chr6 staged outputs"
  pixi_run samtools quickcheck -v "$ILLUMINA_35" "$ILLUMINA_30" "$PACBIO_48" "$PACBIO_30"
  pixi_run bcftools index -n "$PANEL" >/dev/null
  pixi_run bcftools index -n "$TRUTH_VCF" >/dev/null

  echo
  echo "Illumina 30x idxstats:"
  pixi_run samtools idxstats "$ILLUMINA_30" | awk -v c="$CHROM" '$1 == c'
  echo "PacBio HiFi 30x idxstats:"
  pixi_run samtools idxstats "$PACBIO_30" | awk -v c="$CHROM" '$1 == c'

  echo
  echo "Manifest:"
  sed -n '1,5p' "$MANIFEST"
}

main() {
  ensure_pixi
  stage_panel
  stage_map
  stage_truth
  stage_variant_json
  stage_reads
  stage_window_beds
  write_manifest
  validate_outputs
}

main "$@"
