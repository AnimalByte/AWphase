# Phase8E ancestry-weighted candidate v1

## Status

Current best candidate, but not final stable until larger multi-window training and holdout validation.

Stable fallback remains Phase7G.

## What Phase8E adds

Phase8E adds ancestry-weighted donor scoring on top of the Phase7G candidate gate:

- panel sample metadata
- inferred donor population priors
- ancestry-weighted donor support
- ancestry-weighted margin/confidence
- weighted-vs-unweighted donor confidence delta

## Six-window result at threshold 0.50

Phase7G stable:
- total_pred_sites_nonzero: 24607
- total_hamming_errors: 75
- total_switch_errors: 32
- mean_hamming_error_rate: 0.0035849
- mean_switch_error_rate: 0.0022235
- mean_truth_correct_pct: 69.3187

Phase8E ancestry-weighted:
- total_pred_sites_nonzero: 24610
- total_hamming_errors: 74
- total_switch_errors: 30
- mean_hamming_error_rate: 0.0035135
- mean_switch_error_rate: 0.0020448
- mean_truth_correct_pct: 69.3385

## Interpretation

Phase8E is a small but clean improvement over Phase7G:
- slightly higher truth_correct_pct
- fewer hamming errors
- fewer switch errors

Next step:
- Phase8F larger multi-window model
