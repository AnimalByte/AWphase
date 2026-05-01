#!/usr/bin/env bash
set -euo pipefail
cd ~/work/awphase

PYTHON=${PYTHON:-python}
WORKDIR=${WORKDIR:-results/phase4c/chr20_wmec}
mkdir -p "$WORKDIR"

LOCAL_CALLS=${LOCAL_CALLS:-results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv}
OBS_TSV=${OBS_TSV:-results/experiments/chr20_read_frontend/readobs_refaware_v3.obs_debug.tsv}
SCAFFOLD_CALLS=${SCAFFOLD_CALLS:-results/phase3/chr20_internal/scaffold_calls.tsv}
CANDIDATE_SCORES=${CANDIDATE_SCORES:-results/phase4/chr20_mec/candidate_scores.v2.tsv}

PYTHONPATH=python $PYTHON python/awphase_py/build_connectivity_regions_v2.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --candidate-scores-tsv "$CANDIDATE_SCORES" \
  --obs-tsv "$OBS_TSV" \
  --scaffold-calls-tsv "$SCAFFOLD_CALLS" \
  --min-candidate-margin 0.02 \
  --min-obs-confidence 0.025 \
  --min-obs-delta 0.45 \
  --max-sites-per-region 24 \
  --max-region-span-bp 100000 \
  --max-read-link-gap-bp 100000 \
  --out-tsv "$WORKDIR/connectivity_regions.tsv" \
  --out-summary-json "$WORKDIR/connectivity_regions.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/select_active_reads_for_wmec_v1.py \
  --regions-tsv "$WORKDIR/connectivity_regions.tsv" \
  --obs-tsv "$OBS_TSV" \
  --scaffold-calls-tsv "$SCAFFOLD_CALLS" \
  --max-coverage 32 \
  --max-reads-per-region 96 \
  --min-obs-confidence 0.025 \
  --min-obs-delta 0.45 \
  --merge-logodds 6.0 \
  --out-tsv "$WORKDIR/active_reads_wmec.tsv" \
  --out-region-summary-tsv "$WORKDIR/active_reads_wmec.region_summary.tsv" \
  --out-summary-json "$WORKDIR/active_reads_wmec.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_wmec_matrix_v1.py \
  --regions-tsv "$WORKDIR/connectivity_regions.tsv" \
  --active-reads-tsv "$WORKDIR/active_reads_wmec.tsv" \
  --candidate-scores-tsv "$CANDIDATE_SCORES" \
  --out-tsv "$WORKDIR/wmec_matrix.tsv" \
  --out-summary-json "$WORKDIR/wmec_matrix.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/solve_wmec_fpt_v1.py \
  --wmec-matrix-tsv "$WORKDIR/wmec_matrix.tsv" \
  --max-exact-sites 18 \
  --beam-width 256 \
  --min-real-fragments 1 \
  --donor-lambda 0.05 \
  --ranker-lambda 0.04 \
  --anchor-lambda 0.75 \
  --switch-lambda 0.10 \
  --out-tsv "$WORKDIR/wmec_solutions.tsv" \
  --out-region-diagnostics-tsv "$WORKDIR/wmec_region_diagnostics.tsv" \
  --out-summary-json "$WORKDIR/wmec_solutions.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/apply_wmec_regions_v1.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --wmec-solutions-tsv "$WORKDIR/wmec_solutions.tsv" \
  --apply-margin 0.55 \
  --min-real-fragments 1 \
  --min-total-fragments 1 \
  --require-anchor-or-donor \
  --min-abs-donor 0.03 \
  --out-tsv "$WORKDIR/local_calls.phase4c.tsv"

echo "Phase4c wMEC pipeline complete: $WORKDIR"
