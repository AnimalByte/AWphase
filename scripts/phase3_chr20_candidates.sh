#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

SHAPEIT5_COMMON_BIN="${SHAPEIT5_COMMON_BIN:-phase_common}"
SHAPEIT5_RARE_BIN="${SHAPEIT5_RARE_BIN:-phase_rare}"
BEAGLE_JAR="${BEAGLE_JAR:-$HOME/tools/beagle.27Feb25.75f.jar}"

TARGET_BCF="${TARGET_BCF:-data/derived/panel_candidates/chr20/HG002_GRCh38_v5.0q_smvar.chr20.bcf}"
PANEL_BCF="${PANEL_BCF:-data/panels/hgdp_1kg_hg38/subset/panel_chr20_target_sites.bcf}"
REGION="${REGION:-chr20:1-5000000}"
SHAPEIT_MAP="${SHAPEIT_MAP:-data/maps/chr20.b38.gmap.gz}"
THREADS="${THREADS:-8}"
DRY_RUN="${DRY_RUN:-1}"

BEAGLE_MAP=""
beagle_map_candidates=(
  data/maps/beagle_grch38/plink.chr20.GRCh38.map \
  data/maps/beagle_grch38/chr20.map \
  data/maps/beagle_grch38/*chr20*.map
)
for candidate in "${beagle_map_candidates[@]}"; do
  if [[ -e "$candidate" ]]; then
    BEAGLE_MAP="$candidate"
    break
  fi
done

OUTDIR="results/panel_candidates/chr20"
mkdir -p "$OUTDIR"

COMMON_OUT="$OUTDIR/shapeit5.common.bcf"
RARE_OUT="$OUTDIR/shapeit5.rare.bcf"
BEAGLE_PREFIX="$OUTDIR/beagle.chr20"

COMMON_LOG="$OUTDIR/shapeit5.common.log.json"
RARE_LOG="$OUTDIR/shapeit5.rare.log.json"
BEAGLE_LOG="$OUTDIR/beagle.log.json"

CMD_PY=(env PYTHONPATH=python python)

[ -f "$TARGET_BCF" ] || { echo "Missing TARGET_BCF: $TARGET_BCF" >&2; exit 1; }
[ -f "$PANEL_BCF" ] || { echo "Missing PANEL_BCF: $PANEL_BCF" >&2; exit 1; }
[ -f "$BEAGLE_JAR" ] || { echo "Missing BEAGLE_JAR: $BEAGLE_JAR" >&2; exit 1; }

COMMON_ARGS=(
  -m awphase_py.panel_candidates.run_shapeit5_common
  --binary "$SHAPEIT5_COMMON_BIN"
  --input-bcf "$TARGET_BCF"
  --reference-bcf "$PANEL_BCF"
  --region "$REGION"
  --output-bcf "$COMMON_OUT"
  --thread "$THREADS"
  --log "$COMMON_LOG"
)
[ -f "$SHAPEIT_MAP" ] && COMMON_ARGS+=( --map "$SHAPEIT_MAP" )
[ "$DRY_RUN" = "1" ] && COMMON_ARGS+=( --dry-run )

"${CMD_PY[@]}" "${COMMON_ARGS[@]}"

RARE_ARGS=(
  -m awphase_py.panel_candidates.run_shapeit5_rare
  --binary "$SHAPEIT5_RARE_BIN"
  --input-bcf "$TARGET_BCF"
  --scaffold-bcf "$COMMON_OUT"
  --region "$REGION"
  --output-bcf "$RARE_OUT"
  --thread "$THREADS"
  --log "$RARE_LOG"
)
[ -f "$SHAPEIT_MAP" ] && RARE_ARGS+=( --map "$SHAPEIT_MAP" )
[ "$DRY_RUN" = "1" ] && RARE_ARGS+=( --dry-run )

"${CMD_PY[@]}" "${RARE_ARGS[@]}"

BEAGLE_ARGS=(
  -m awphase_py.panel_candidates.run_beagle_phase
  --jar "$BEAGLE_JAR"
  --gt "$TARGET_BCF"
  --ref "$PANEL_BCF"
  --out-prefix "$BEAGLE_PREFIX"
  --threads "$THREADS"
  --log "$BEAGLE_LOG"
)
[ -n "$BEAGLE_MAP" ] && BEAGLE_ARGS+=( --map "$BEAGLE_MAP" )
[ "$DRY_RUN" = "1" ] && BEAGLE_ARGS+=( --dry-run )

"${CMD_PY[@]}" "${BEAGLE_ARGS[@]}"

echo
echo "Phase 3A chr20 candidate-generation run complete."
echo "Target BCF:  $TARGET_BCF"
echo "Panel BCF:   $PANEL_BCF"
echo "Shapeit map: ${SHAPEIT_MAP:-NONE}"
echo "Beagle map:  ${BEAGLE_MAP:-NONE}"
echo "Logs:"
echo "  $COMMON_LOG"
echo "  $RARE_LOG"
echo "  $BEAGLE_LOG"
