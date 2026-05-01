#!/usr/bin/env bash
set -euo pipefail
cd ~/work/awphase

WORKDIR=${WORKDIR:-results/phase4d/chr20_orient}
mkdir -p "$WORKDIR"

PYTHONPATH=python python python/awphase_py/orient_wmec_regions_v1.py \
  --regions-tsv results/phase4c/chr20_wmec/connectivity_regions.tsv \
  --wmec-solutions-tsv results/phase4c/chr20_wmec/wmec_solutions.tsv \
  --local-calls-tsv results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv \
  --obs-tsv results/experiments/chr20_read_frontend/readobs_refaware_v3.obs_debug.tsv \
  --out-oriented-tsv "$WORKDIR/region_resolutions.oriented.tsv" \
  --out-orientations-tsv "$WORKDIR/region_orientations.tsv" \
  --out-summary-json "$WORKDIR/region_orientations.summary.json"

PYTHONPATH=python python python/awphase_py/apply_region_resolved_calls_v1.py \
  --local-calls-tsv results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv \
  --region-resolutions-tsv "$WORKDIR/region_resolutions.oriented.tsv" \
  --apply-margin 0.75 \
  --min-query-margin 0.15 \
  --require-anchor-or-donor \
  --out-tsv "$WORKDIR/local_calls.phase4d.tsv"

echo "Phase 4d orientation pass complete: $WORKDIR"
