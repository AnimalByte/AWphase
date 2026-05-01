#!/usr/bin/env bash
set -euo pipefail
cd ~/work/awphase

PYTHON=${PYTHON:-python}
CHROM=${CHROM:-chr20}

LOCAL_CALLS=${LOCAL_CALLS:-results/experiments/chr20_refaware_v3_selected_on/local_calls.tsv}
SELECTED_READS=${SELECTED_READS:-results/experiments/chr20_read_frontend/selected_reads.tsv}
HAPLOTAG=${HAPLOTAG:-results/experiments/chr20_read_frontend/haplotag_lite.tsv}
SUPERREADS=${SUPERREADS:-results/experiments/chr20_read_frontend/superreads.tsv}
OBS_TSV=${OBS_TSV:-results/experiments/chr20_read_frontend/readobs_refaware_v3.obs_debug.tsv}
SITE_DEBUG=${SITE_DEBUG:-results/experiments/chr20_read_frontend/readobs_refaware_v3.site_debug.tsv}

SITE_PRIORS=${SITE_PRIORS:-data/derived/hg002_chr20_site_priors_ancestry_on.json}
PANEL_CANDIDATES=${PANEL_CANDIDATES:-data/derived/panel_chr20_candidates.json}
TARGET_SITES=${TARGET_SITES:-data/panels/hgdp_1kg_hg38/subset/hg002_chr20_target_sites.tsv}
VARIANT_JSON=${VARIANT_JSON:-data/derived/hg002_chr20_variants.real.json}
TRUTH_VCF=${TRUTH_VCF:-data/truth/hg002_phase_chr20/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.chr20.vcf.gz}
FASTA=${FASTA:-references/grch38_noalt_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}

WORKDIR=${WORKDIR:-results/phase3/chr20_internal}
mkdir -p "$WORKDIR"

PYTHONPATH=python $PYTHON python/awphase_py/build_scaffold_common_v2.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --selected-reads-tsv "$SELECTED_READS" \
  --haplotag-tsv "$HAPLOTAG" \
  --superreads-tsv "$SUPERREADS" \
  --obs-tsv "$OBS_TSV" \
  --site-debug-tsv "$SITE_DEBUG" \
  --site-priors-json "$SITE_PRIORS" \
  --chrom "$CHROM" \
  --out-calls-tsv "$WORKDIR/scaffold_calls.tsv" \
  --out-positions-tsv "$WORKDIR/scaffold_positions.tsv" \
  --out-summary-json "$WORKDIR/scaffold_summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_candidate_queries_v1.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --scaffold-calls-tsv "$WORKDIR/scaffold_calls.tsv" \
  --out-tsv "$WORKDIR/candidate_queries.tsv" \
  --out-summary-json "$WORKDIR/candidate_queries.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_donor_prior_v1.py \
  --candidate-queries-tsv "$WORKDIR/candidate_queries.tsv" \
  --panel-candidates-json "$PANEL_CANDIDATES" \
  --target-sites-tsv "$TARGET_SITES" \
  --scaffold-calls-tsv "$WORKDIR/scaffold_calls.tsv" \
  --window-sites 12 \
  --min-anchor-count 2 \
  --top-groups 2 \
  --out-tsv "$WORKDIR/donor_prior.tsv" \
  --out-summary-json "$WORKDIR/donor_prior.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_reference_context_features_v1.py \
  --variant-json "$VARIANT_JSON" \
  --fasta "$FASTA" \
  --chrom "$CHROM" \
  --flank 12 \
  --out-tsv "$WORKDIR/reference_context.tsv" \
  --out-summary-json "$WORKDIR/reference_context.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_candidate_features_v1.py \
  --candidate-queries-tsv "$WORKDIR/candidate_queries.tsv" \
  --scaffold-calls-tsv "$WORKDIR/scaffold_calls.tsv" \
  --donor-prior-tsv "$WORKDIR/donor_prior.tsv" \
  --reference-context-tsv "$WORKDIR/reference_context.tsv" \
  --out-tsv "$WORKDIR/candidate_features.tsv" \
  --out-summary-json "$WORKDIR/candidate_features.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/build_candidate_labels_from_truth_v1.py \
  --candidate-features-tsv "$WORKDIR/candidate_features.tsv" \
  --variant-json "$VARIANT_JSON" \
  --truth-vcf "$TRUTH_VCF" \
  --out-tsv "$WORKDIR/candidate_features.labeled.tsv" \
  --out-summary-json "$WORKDIR/candidate_features.labeled.summary.json"

PYTHONPATH=python $PYTHON python/awphase_py/train_candidate_ranker_v1.py \
  --training-tsv "$WORKDIR/candidate_features.labeled.tsv" \
  --out-model-json "$WORKDIR/candidate_ranker_v1.json" \
  --out-metrics-json "$WORKDIR/candidate_ranker_v1.metrics.json" \
  --out-feature-importance-tsv "$WORKDIR/candidate_ranker_v1.feature_importance.tsv"

PYTHONPATH=python $PYTHON python/awphase_py/score_candidate_ranker_v1.py \
  --inference-tsv "$WORKDIR/candidate_features.tsv" \
  --model-json "$WORKDIR/candidate_ranker_v1.json" \
  --out-tsv "$WORKDIR/candidate_scores.tsv"

PYTHONPATH=python $PYTHON python/awphase_py/apply_ranked_candidates_v1.py \
  --local-calls-tsv "$LOCAL_CALLS" \
  --scored-candidates-tsv "$WORKDIR/candidate_scores.tsv" \
  --out-tsv "$WORKDIR/local_calls.ranked.tsv" \
  --min-best-score 0.60 \
  --min-margin 0.05 \
  --min-modify-score 0.72 \
  --min-modify-margin 0.15 \
  --anchor-agreement-min 0.0 \
  --strong-isolated-score 0.85 \
  --strong-isolated-margin 0.25 \
  --max-gap-bp 5000

echo "Phase 3 internal pipeline complete: $WORKDIR"
