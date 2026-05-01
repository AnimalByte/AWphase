# Phase8F ancestry-weighted multi-window candidate v1

## Status

Current best learned candidate.

## Training set

17 HG002 windows across chr1, chr20, and chr22.

## Main result on same 17 windows

Phase7A hand panel fill:
- total_pred_sites_nonzero: 64837
- total_hamming_errors: 300
- total_switch_errors: 263
- mean_truth_correct_pct: 70.8709

Phase8F threshold 0.50:
- total_pred_sites_nonzero: 64413
- total_hamming_errors: 250
- total_switch_errors: 207
- mean_truth_correct_pct: 70.6001

Interpretation:
- Phase8F t0.50 keeps nearly the same recall as Phase7A
- reduces hamming errors by about 16.7%
- reduces switch errors by about 21.3%

Recommended modes:
- t0.50: balanced/default candidate
- t0.80: precision mode
- t0.90: high-precision mode
