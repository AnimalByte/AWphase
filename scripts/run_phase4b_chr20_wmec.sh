#!/usr/bin/env bash
set -euo pipefail
cd ~/work/awphase

PYTHON=${PYTHON:-python}
WORKDIR=${WORKDIR:-results/phase4b/chr20_wmec}
mkdir -p "$WORKDIR"

# Reuse existing chr20 artifacts
LOCAL_CALLS=${LOCAL_CALLS:-results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv}
OBS_TSV=${OBS_TSV:-results/experiments/chr20_read_frontend/readobs_refaware_v3.obs_debug.tsv}
SCAFFOLD_CALLS=${SCAFFOLD_CALLS:-results/phase3/chr20_internal/scaffold_calls.tsv}
CANDIDATE_SCORES=${CANDIDATE_SCORES:-results/phase4/chr20_mec/candidate_scores.v2.tsv}

PYTHONPATH=python $PYTHON python/awphase_py/build_connectivity_regions_v2.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --candidate-scores-tsv "$CANDIDATE_SCORES" \
  --obs-tsv "$OBS_TSV" \
  --scaffold-calls-tsv "$SCAFFOLD_CALLS" \
  --out-tsv "$WORKDIR/connectivity_regions.tsv" \
  --out-summary-json "$WORKDIR/connectivity_regions.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/select_active_reads_for_wmec_v1.py \
  --regions-tsv "$WORKDIR/connectivity_regions.tsv" \
  --obs-tsv "$OBS_TSV" \
  --scaffold-calls-tsv "$SCAFFOLD_CALLS" \
  --out-tsv "$WORKDIR/active_reads_wmec.tsv" \
  --out-summary-json "$WORKDIR/active_reads_wmec.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_wmec_matrix_v1.py \
  --regions-tsv "$WORKDIR/connectivity_regions.tsv" \
  --active-reads-tsv "$WORKDIR/active_reads_wmec.tsv" \
  --candidate-scores-tsv "$CANDIDATE_SCORES" \
  --out-tsv "$WORKDIR/wmec_matrix.tsv" \
  --out-summary-json "$WORKDIR/wmec_matrix.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/solve_wmec_fpt_v1.py \
  --wmec-matrix-tsv "$WORKDIR/wmec_matrix.tsv" \
  --out-tsv "$WORKDIR/wmec_solutions.tsv" \
  --out-summary-json "$WORKDIR/wmec_solutions.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/apply_wmec_regions_v1.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --wmec-solutions-tsv "$WORKDIR/wmec_solutions.tsv" \
  --apply-margin 1.0 \
  --require-anchor-or-donor \
  --out-tsv "$WORKDIR/local_calls.phase4b.tsv"

echo "Phase4b wMEC pipeline complete: $WORKDIR"
