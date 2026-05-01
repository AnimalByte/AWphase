#!/usr/bin/env bash
set -euo pipefail
cd ~/work/awphase

PYTHON=${PYTHON:-python}
CHROM=${CHROM:-chr20}
LOCAL_CALLS=${LOCAL_CALLS:-results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv}
SELECTED_READS=${SELECTED_READS:-results/experiments/chr20_read_frontend/selected_reads.tsv}
HAPLOTAG=${HAPLOTAG:-results/experiments/chr20_read_frontend/haplotag_lite.tsv}
OBS_TSV=${OBS_TSV:-results/experiments/chr20_read_frontend/readobs_refaware_v3.obs_debug.tsv}
SITE_PRIORS=${SITE_PRIORS:-data/derived/hg002_chr20_site_priors_ancestry_on.json}
PANEL_CANDIDATES=${PANEL_CANDIDATES:-data/derived/panel_chr20_candidates.json}
TARGET_SITES=${TARGET_SITES:-data/panels/hgdp_1kg_hg38/subset/hg002_chr20_target_sites.tsv}
VARIANT_JSON=${VARIANT_JSON:-data/derived/hg002_chr20_variants.real.json}
TRUTH_VCF=${TRUTH_VCF:-data/truth/hg002_phase_chr20/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.chr20.vcf.gz}
BASE=${BASE:-results/phase3/chr20_internal}
WORKDIR=${WORKDIR:-results/phase4/chr20_mec}
mkdir -p "$WORKDIR"

if [[ ! -f "$BASE/scaffold_calls.tsv" || ! -f "$BASE/candidate_queries.tsv" || ! -f "$BASE/candidate_features.tsv" ]]; then
  bash scripts/run_phase3_chr20_internal.sh
fi

PYTHONPATH=python $PYTHON python/awphase_py/build_scaffold_as_reads_v1.py   --scaffold-calls-tsv "$BASE/scaffold_calls.tsv"   --out-tsv "$WORKDIR/scaffold_as_reads.tsv"   --out-summary-json "$WORKDIR/scaffold_as_reads.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_local_donor_prior_v2.py   --candidate-queries-tsv "$BASE/candidate_queries.tsv"   --panel-candidates-json "$PANEL_CANDIDATES"   --target-sites-tsv "$TARGET_SITES"   --scaffold-calls-tsv "$BASE/scaffold_calls.tsv"   --out-tsv "$WORKDIR/donor_prior_v2.tsv"   --out-summary-json "$WORKDIR/donor_prior_v2.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_realign_like_features_v1.py   --candidate-queries-tsv "$BASE/candidate_queries.tsv"   --obs-tsv "$OBS_TSV"   --haplotag-tsv "$HAPLOTAG"   --out-tsv "$WORKDIR/realign_like_features.tsv"   --out-summary-json "$WORKDIR/realign_like_features.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_candidate_features_v2.py   --base-features-tsv "$BASE/candidate_features.tsv"   --donor-prior-v2-tsv "$WORKDIR/donor_prior_v2.tsv"   --realign-features-tsv "$WORKDIR/realign_like_features.tsv"   --scaffold-calls-tsv "$BASE/scaffold_calls.tsv"   --out-tsv "$WORKDIR/candidate_features.v2.tsv"   --out-summary-json "$WORKDIR/candidate_features.v2.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_candidate_labels_from_truth_v1.py   --candidate-features-tsv "$WORKDIR/candidate_features.v2.tsv"   --variant-json "$VARIANT_JSON"   --truth-vcf "$TRUTH_VCF"   --out-tsv "$WORKDIR/candidate_features.v2.labeled.tsv"   --out-summary-json "$WORKDIR/candidate_features.v2.labeled.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/train_candidate_ranker_v2.py   --training-tsv "$WORKDIR/candidate_features.v2.labeled.tsv"   --out-model-json "$WORKDIR/candidate_ranker_v2.json"   --out-metrics-json "$WORKDIR/candidate_ranker_v2.metrics.json"   --out-feature-importance-tsv "$WORKDIR/candidate_ranker_v2.feature_importance.tsv"

PYTHONPATH=python $PYTHON python/awphase_py/score_candidate_ranker_v2.py   --inference-tsv "$WORKDIR/candidate_features.v2.tsv"   --model-json "$WORKDIR/candidate_ranker_v2.json"   --out-tsv "$WORKDIR/candidate_scores.v2.tsv"

PYTHONPATH=python $PYTHON python/awphase_py/build_region_queries_v1.py   --candidate-scores-tsv "$WORKDIR/candidate_scores.v2.tsv"   --scaffold-calls-tsv "$BASE/scaffold_calls.tsv"   --out-tsv "$WORKDIR/region_queries.tsv"   --out-summary-json "$WORKDIR/region_queries.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/resolve_regions_mec_v1.py   --region-queries-tsv "$WORKDIR/region_queries.tsv"   --candidate-scores-tsv "$WORKDIR/candidate_scores.v2.tsv"   --obs-tsv "$OBS_TSV"   --scaffold-as-reads-tsv "$WORKDIR/scaffold_as_reads.tsv"   --out-tsv "$WORKDIR/region_resolutions.tsv"   --out-summary-json "$WORKDIR/region_resolutions.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/apply_region_resolved_calls_v1.py   --local-calls-tsv "$LOCAL_CALLS"   --region-resolutions-tsv "$WORKDIR/region_resolutions.tsv"   --apply-margin 0.75   --min-query-margin 0.15   --require-anchor-or-donor   --out-tsv "$WORKDIR/local_calls.phase4.tsv"

echo "Phase4 MEC-like pipeline complete: $WORKDIR"
