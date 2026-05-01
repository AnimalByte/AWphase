#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-8}"
CURL_RETRIES="${CURL_RETRIES:-5}"
PIXI_BIN="${PIXI_BIN:-$HOME/.pixi/bin/pixi}"
CHROMS="${CHROMS:-chr1 chr20 chr22 chr6}"
INCLUDE_REFERENCE="${INCLUDE_REFERENCE:-1}"
read -r -a CHROM_ARRAY <<< "$CHROMS"

GNOMAD_PANEL_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
GNOMAD_META_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv"
BEAGLE_MAP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"
GIAB_REF_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38"

REF_GZ="references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
REF_FASTA="references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
MAP_ZIP="data/maps/beagle_grch38/plink.GRCh38.map.zip"
MAP_DIR="data/maps/beagle_grch38/no_chr_in_chrom_field"
PANEL_DIR="data/panels/hgdp_1kg_hg38/full"
META_DIR="data/panels/hgdp_1kg_hg38/meta"

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

stage_panels() {
  log "Stage HGDP + 1KGP SHAPEIT5 reference panels"
  mkdir -p "$PANEL_DIR" "$META_DIR"
  for chrom in "${CHROM_ARRAY[@]}"; do
    local base="hgdp1kgp_${chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
    local panel="${PANEL_DIR}/${base}"
    download_url "${GNOMAD_PANEL_BASE}/${base}" "$panel"
    download_url "${GNOMAD_PANEL_BASE}/${base}.csi" "${panel}.csi"
    pixi_run bcftools index -n "$panel" >/dev/null
  done
  download_url "$GNOMAD_META_URL" "${META_DIR}/gnomad_meta_updated.tsv"
}

stage_maps() {
  log "Stage Beagle GRCh38 genetic maps"
  download_url "$BEAGLE_MAP_URL" "$MAP_ZIP"
  mkdir -p "$MAP_DIR"
  pixi_run python - "$MAP_ZIP" "$MAP_DIR" "${CHROM_ARRAY[@]}" <<'PY'
import sys
import zipfile
from pathlib import Path

zip_path = sys.argv[1]
map_dir = Path(sys.argv[2])
chroms = sys.argv[3:]

with zipfile.ZipFile(zip_path) as zf:
    for chrom in chroms:
        member = f"no_chr_in_chrom_field/plink.{chrom}.GRCh38.map"
        out = map_dir / f"plink.{chrom}.GRCh38.map"
        if out.exists() and out.stat().st_size:
            print(f"already present: {out}")
            continue
        data = zf.read(member)
        out.write_bytes(data)
        print({"out": str(out), "bytes": len(data)})
PY
}

stage_reference() {
  if [[ "$INCLUDE_REFERENCE" != "1" ]]; then
    echo "SKIP reference download because INCLUDE_REFERENCE=$INCLUDE_REFERENCE"
    return
  fi

  log "Stage GRCh38 no-alt reference"
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

main() {
  ensure_pixi
  stage_panels
  stage_maps
  stage_reference

  echo
  echo "Runtime data staged for chromosomes: $CHROMS"
  echo "Panels: $PANEL_DIR"
  echo "Maps:   $MAP_DIR"
  if [[ "$INCLUDE_REFERENCE" == "1" ]]; then
    echo "Ref:    $REF_FASTA"
  fi
}

main "$@"
