#!/usr/bin/env bash
set -euo pipefail

BASECFG="configs/benchmark_chr20.tuned.yaml"
OUTROOT="results/sweeps/bridge_thresholds_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTROOT/configs"

SUMMARY="$OUTROOT/summary.tsv"
echo -e "min_confidence\tmin_informative_sites\tnumber_phased\tnumber_abstained\tjoined_boundaries\tflipped_boundaries\tabstained_boundaries\tn_exact_overlap_sites_phased\thamming_error_rate\tswitch_error_rate\tphased_site_accuracy_pct\ttruth_correct_pct\traw_block_n50_bp\tcfg_path\tmetrics_json" > "$SUMMARY"

for conf in 0.15 0.25 0.35 0.45 0.55 0.65 0.75; do
  for info in 3 5 7 10; do
    tag="conf_${conf}_info_${info}"
    cfg="$OUTROOT/configs/${tag}.yaml"
    run_dir="$OUTROOT/$tag"
    mkdir -p "$run_dir"

    cp "$BASECFG" "$cfg"

    sed -i -E "s/^([[:space:]]*min_confidence:).*/\1 ${conf}/" "$cfg"
    sed -i -E "s/^([[:space:]]*min_informative_sites:).*/\1 ${info}/" "$cfg"

    echo
    echo "== running $tag =="

    ./target/debug/awphase-cli --config "$cfg" benchmark > "$run_dir/benchmark.stdout.json"

    cp results/runs/HG002/ancestry_on/benchmark_summary.json        "$run_dir/benchmark_summary.json"
    cp results/runs/HG002/ancestry_on/boundary_decisions.json       "$run_dir/boundary_decisions.json"
    cp results/runs/HG002/ancestry_on/derived_boundary_scores.json  "$run_dir/derived_boundary_scores.json"
    cp results/runs/HG002/ancestry_on/derived_block_summaries.json  "$run_dir/derived_block_summaries.json"
    cp results/runs/HG002/ancestry_on/final_stitched_calls.tsv      "$run_dir/final_stitched_calls.tsv"

    python python/awphase_py/evaluate_phase_truth.py \
      --pred-tsv "$run_dir/final_stitched_calls.tsv" \
      --variant-json data/derived/hg002_chr20_variants.real.json \
      --truth-vcf data/truth/hg002_phase_chr20/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.chr20.vcf.gz \
      --truth-bed data/truth/hg002_phase_chr20/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed \
      --phase-column stitched_phase_state \
      --out-prefix "$run_dir/truth_eval_stitched" \
      > "$run_dir/eval.stdout.json"

    python - "$conf" "$info" "$run_dir/benchmark_summary.json" "$run_dir/truth_eval_stitched.metrics.json" "$cfg" "$SUMMARY" <<'PY'
import sys, json
conf, info, bench_path, metrics_path, cfg_path, summary_path = sys.argv[1:]

with open(bench_path) as fh:
    b = json.load(fh)
with open(metrics_path) as fh:
    m = json.load(fh)

row = [
    conf,
    info,
    str(b.get("number_phased", "")),
    str(b.get("number_abstained", "")),
    str(b.get("joined_boundaries", "")),
    str(b.get("flipped_boundaries", "")),
    str(b.get("abstained_boundaries", "")),
    str(m.get("n_exact_overlap_sites_phased", "")),
    str(m.get("hamming_error_rate", "")),
    str(m.get("switch_error_rate", "")),
    str(m.get("phased_site_accuracy_pct", "")),
    str(m.get("truth_correct_pct", "")),
    str(m.get("raw_block_n50_bp", "")),
    cfg_path,
    metrics_path,
]

with open(summary_path, "a") as out:
    out.write("\t".join(row) + "\n")

print("\t".join(row[:13]))
PY

  done
done

echo
echo "Sweep complete."
echo "Summary: $SUMMARY"
