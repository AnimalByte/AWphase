#!/usr/bin/env bash
set -euo pipefail

cd ~/work/awphase

OUT="results/comparisons/chr22_awphase_vs_real_whatshap"
mkdir -p "$OUT"

SUMMARY="$OUT/summary.tsv"

python - <<'PY'
import json
from pathlib import Path

branches = {
    # AWPhase internal
    "awphase_tagged": "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d1000/eval.metrics.json",
    "awphase_readbridge_superreads_d500": "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d500/eval.metrics.json",
    "awphase_readbridge_superreads_d1000": "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d1000/eval.metrics.json",
    "awphase_readbridge_superreads_d2500": "results/phase5b/chr22_20_25mb_readbridge_superreads_safe/graft_d2500/eval.metrics.json",
    "awphase_readbridge_selected_d5000": "results/phase5b/chr22_20_25mb_readbridge_selected_safe/graft_d5000/eval.metrics.json",

    # Real WhatsHap generated here
    "whatshap_illumina_30x_chr22": "baselines/chr22_20_25mb/illumina_30x/truth_eval.metrics.json",
    "whatshap_pacbio_30x_chr22": "baselines/chr22_20_25mb/pacbio_30x/truth_eval.metrics.json",
}

# Correct tagged baseline if available from earlier direct evaluation.
tagged_candidates = [
    "results/phase5b/chr22_20_25mb_distance_gated/tagged.eval.metrics.json",
    "results/phase5b/chr22_20_25mb_safe/tagged_baseline.eval.metrics.json",
]
for p in tagged_candidates:
    if Path(p).exists():
        branches["awphase_tagged_baseline"] = p
        break

keys = [
    "n_pred_sites_nonzero",
    "n_exact_overlap_sites_phased",
    "hamming_denominator",
    "hamming_error_rate",
    "switch_denominator",
    "switch_error_rate",
    "phased_site_accuracy_pct",
    "truth_correct_pct",
    "raw_block_n50_bp",
    "max_block_span_bp",
    "median_block_span_bp",
]

rows = []
for name, path in branches.items():
    p = Path(path)
    if not p.exists():
        print(f"missing: {name} -> {path}")
        continue
    m = json.load(open(p))
    rows.append((name, path, m))

out = Path("results/comparisons/chr22_awphase_vs_real_whatshap/summary.tsv")
with out.open("w") as fh:
    fh.write("branch\tpath\t" + "\t".join(keys) + "\n")
    for name, path, m in rows:
        fh.write(name + "\t" + path + "\t" + "\t".join(str(m.get(k, "")) for k in keys) + "\n")

print("Wrote", out)
PY

echo
echo "===== all branches ====="
column -t "$SUMMARY"

echo
echo "===== best hamming ====="
column -t "$SUMMARY" | sort -k6,6n | sed -n '1,30p'

echo
echo "===== best switch ====="
column -t "$SUMMARY" | sort -k8,8n | sed -n '1,30p'
