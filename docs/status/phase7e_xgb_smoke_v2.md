# Phase7E XGBoost acceptance smoke model v2

Historical note: this smoke-model note is not the current repository default.
The current default is tracked in `docs/current_default.md`.

## Summary

Phase7E trains an XGBoost accept/reject model for Phase7A panel-fill candidates.

Training data:
- rows: 1016
- positives: 893
- negatives: 123
- physical groups: chr20_30_35mb, chr22_20_25mb

Grouped CV collapses chr20 30x and chr20 35x into one physical group.

## Main result

XGBoost provides a useful precision/recall gate.

Recommended operating modes:

- Phase7E-recall: threshold 0.50
- Phase7E-balanced: threshold 0.70
- Phase7E-precision: threshold 0.80

Threshold 0.70 was the best default for this smoke-model milestone.

## Gated evaluation

chr22:
- t0.70: truth_correct_pct 66.21, hamming 4, switch 4
- t0.80: truth_correct_pct 64.99, hamming 1, switch 1

chr20 30x:
- t0.70: truth_correct_pct 58.27, hamming 6, switch 6
- t0.80: truth_correct_pct 57.33, hamming 5, switch 5

chr20 35x:
- t0.70: truth_correct_pct 58.75, hamming 6, switch 6
- t0.80: truth_correct_pct 57.88, hamming 5, switch 5

## Feature importance

Strongest signal came from panel_margin, panel_confidence, best_vs_second_margin, panel_conflict, and support_conflict_ratio. candidate_accepted_by_rule was not dominant, so the model is not simply copying the hand threshold.

## Interpretation

This is a smoke model, not final. It has only two independent physical groups. More windows are needed before finalizing the Phase7E acceptance model.
