# Phase7E XGBoost v3

Historical note: this document describes the Phase7E recommendation at the time
it was written. The current repository default is tracked in
`docs/current_default.md`.

## Summary

Phase7E v3 is a learned XGBoost accept/reject gate for Phase7A panel-fill candidates.

Phase7A uses hand thresholds to accept panel-LD fills. Phase7E learns a smoother acceptance surface from labeled candidate fills using truth-oriented labels.

## Training data

- Rows: 1,733
- Positives: 1,527
- Negatives: 206
- Physical groups: 4
  - chr20_25_30mb
  - chr20_30_35mb
  - chr20_45_50mb
  - chr22_20_25mb

The chr20 30x and 35x versions are collapsed into one physical group for cross-validation.

## Cross-validation

- Grouped CV folds: 4
- CV ROC-AUC: 0.7677
- CV average precision: 0.9615
- CV accuracy at threshold 0.5: 0.7149

## Feature importance

The model is not copying the original hand threshold. `candidate_accepted_by_rule` has zero importance.

Most useful features:

- panel_margin
- panel_confidence
- support_conflict_ratio
- panel_support
- panel_conflict
- best_vs_second_margin
- n_block_candidates

## Recommended Phase7E modes

- Phase7E-recall: XGB threshold 0.50
- Phase7E-balanced: XGB threshold 0.70
- Phase7E-precision: XGB threshold 0.80

Historical default recommendation for this milestone: Phase7E-balanced, threshold 0.70.

## Current interpretation

AWPhase Phase6C is a high-precision, lower-recall phaser.

Phase7A improves recall by adding panel-LD fills, but introduces additional hamming/switch errors.

Phase7E v3 reduces many Phase7A errors while preserving much of the recall gain.

## Baselines

WhatsHap Illumina 30x generally has higher recall than AWPhase but more switch/hamming errors in these chr20 windows.

WhatsHap PacBio 30x has very high truth-correct recall and long blocks, but the current PacBio baseline needs more validation. Orientation-aware diagnostics show residual wrong sites remain, especially in chr20_45_50mb, where block `whatshap_46493090` has a large residual error cluster.

## Next steps

1. Add chr1 and additional chr22 30x Illumina windows.
2. Add WhatsHap Illumina baselines for those windows.
3. Add PacBio baselines where PacBio data exists or can be sliced.
4. Retrain Phase7E v4 with more chromosome diversity.
5. Investigate PacBio chr20_45_50mb block-level residual errors.
