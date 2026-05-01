# Phase 6A chr22 20-25Mb wMEC status

## Summary

Phase 6A replaces the Phase 5 graft-centered layer with direct read-fragment weighted MEC solving.

## Inputs

- obs/debug: results/experiments/chr22_refaware_v3_on/readobs_refaware_v3.obs_debug.tsv
- variants: data/derived/hg002_chr22_variants.real.json
- truth: HG002 v5.0q chr22 20-25Mb

## Phase 6A metrics

- n_pred_sites_nonzero: 3068
- hamming_error_rate: 0.0009302325581395349
- switch_error_rate: 0.0012461059190031153
- phased_site_accuracy_pct: 99.90697674418605
- truth_correct_pct: 60.643704121964994
- raw_block_n50_bp: 1770
- max_block_span_bp: 146784
- median_block_span_bp: 522

## Comparison

Phase 6A is dramatically more accurate than Phase 5b and locally more accurate than WhatsHap Illumina 30x on this window, but it is much more fragmented.

WhatsHap Illumina 30x:
- hamming_error_rate: 0.008959044368600682
- switch_error_rate: 0.009313154831199068
- truth_correct_pct: 65.5844155844156
- raw_block_n50_bp: 1336301

## Diagnostic conclusion

Phase 5b failed due to mixed/noisy within-block patterns.
Phase 6A fixes local phase coherence.

Block-mixing vs WhatsHap:
- blocks_evaluated: 689
- residual_error_rate_after_best_block_orientation: 0.0009778357235984355
- adjacent_pair_mismatch_rate: 0.0008406893652795292
- mixed_blocks: 2

## Next direction

Phase 6B should stitch small high-confidence MEC components into longer phase sets while preserving Phase 6A local correctness.
