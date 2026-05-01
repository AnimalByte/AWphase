#!/usr/bin/env bash
set -u
set +e
set +o pipefail

cd ~/work/awphase || exit 2

mkdir -p results/phase9/logs

COMBINE_LOG="results/phase9/logs/phase9d_combine_wide_noleak.crashsafe.log"
TRAIN_LOG="results/phase9/logs/phase9d_train_wide_xgb_noleak.crashsafe.log"

TRAIN_TSV="results/phase9/training.candidates.phase9d.block_bridges.wide.noleak.tsv"
PRED_TSV="results/phase9/phase9d_wide_bridge_xgb.noleak.oof_predictions.tsv"

echo
echo "===== compile trainer ====="
python -m py_compile python/awphase_py/phase9/train_phase9d_wide_bridge_xgb_noleak_safe.py
rc=$?
if [[ "$rc" != "0" ]]; then
  echo "compile failed rc=$rc"
  exit "$rc"
fi

echo
echo "===== combine wide candidates ====="
python scripts/phase9/combine_phase9d_wide_noleak_safe.py > "$COMBINE_LOG" 2>&1
rc=$?

echo "combine rc=$rc"
tail -n 160 "$COMBINE_LOG"

if [[ "$rc" != "0" ]]; then
  echo
  echo "COMBINE FAILED. You probably need to generate Phase9D wide candidates first."
  exit "$rc"
fi

if [[ ! -s "$TRAIN_TSV" ]]; then
  echo "COMBINE FAILED: expected training table missing: $TRAIN_TSV"
  exit 2
fi

echo
echo "===== training table quick check ====="
python - <<PY
import pandas as pd
from pathlib import Path
p = Path("$TRAIN_TSV")
print("file:", p)
print("size_mb:", round(p.stat().st_size / 1024 / 1024, 2))
df = pd.read_csv(p, sep="\\t", nrows=5)
print("n_cols:", len(df.columns))
print("columns:")
for c in df.columns:
    print("  " + c)
PY

echo
echo "===== train Phase9D crash-safe model ====="

PYTHONPATH=python python -u python/awphase_py/phase9/train_phase9d_wide_bridge_xgb_noleak_safe.py \
  --training-tsv "$TRAIN_TSV" \
  --out-model-json results/phase9/phase9d_wide_bridge_xgb.noleak.json \
  --out-metrics-json results/phase9/phase9d_wide_bridge_xgb.noleak.metrics.json \
  --out-feature-importance-tsv results/phase9/phase9d_wide_bridge_xgb.noleak.feature_importance.tsv \
  --out-predictions-tsv "$PRED_TSV" \
  --out-feature-report-json results/phase9/phase9d_wide_bridge_xgb.noleak.feature_report.json \
  --threads 6 \
  --n-splits 5 \
  --n-estimators 250 \
  --max-depth 3 \
  --learning-rate 0.04 \
  > "$TRAIN_LOG" 2>&1

rc=$?

echo
echo "===== train rc=$rc ====="
tail -n 260 "$TRAIN_LOG"

if [[ "$rc" != "0" ]]; then
  echo
  echo "TRAINING FAILED. See full log:"
  echo "$TRAIN_LOG"
  exit "$rc"
fi

if [[ ! -s "$PRED_TSV" ]]; then
  echo "TRAINING CLAIMED SUCCESS BUT PREDICTIONS ARE MISSING:"
  echo "$PRED_TSV"
  exit 3
fi

echo
echo "===== output files ====="
ls -lh \
  results/phase9/phase9d_wide_bridge_xgb.noleak.json \
  results/phase9/phase9d_wide_bridge_xgb.noleak.metrics.json \
  results/phase9/phase9d_wide_bridge_xgb.noleak.feature_importance.tsv \
  results/phase9/phase9d_wide_bridge_xgb.noleak.feature_report.json \
  "$PRED_TSV"

echo
echo "===== metrics head ====="
python - <<'PY'
import json
m = json.load(open("results/phase9/phase9d_wide_bridge_xgb.noleak.metrics.json"))
for k in [
    "rows",
    "positives",
    "negatives",
    "positive_rate",
    "selected_feature_count",
    "cv_roc_auc",
    "cv_average_precision",
    "cv_accuracy_0.5",
    "cv_precision_0.5",
    "cv_recall_0.5",
    "cv_f1_0.5",
]:
    print(f"{k}: {m.get(k)}")
PY

echo
echo "===== feature importance ====="
column -t -s $'\t' results/phase9/phase9d_wide_bridge_xgb.noleak.feature_importance.tsv | head -n 80
