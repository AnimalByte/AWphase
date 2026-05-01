# Phase 5b: Tagged-spine nearest-block graft

## Status

Current best chr20 AWPhase branch.

## Method

Uses the `tagged` branch as the continuity spine, then grafts Phase 5/Phase 4e filled calls into the nearest tagged phase block. This fixes the earlier issue where filled calls increased `n_pred_sites_nonzero` but did not enter the hamming/switch denominator.

## chr20 v5q result

- n_pred_sites_nonzero: 4734
- hamming_error_rate: 0.15214592274678113
- switch_error_rate: 0.16094541361845807
- phased_site_accuracy_pct: 84.78540772532187
- truth_correct_pct: 6.79519812878371
- raw_block_n50_bp: 6279

## Interpretation

This branch beats the tagged baseline on hamming, switch, phased-site accuracy, and truth_correct_pct on chr20 v5q. Next step is cross-chrom validation on chr1 and chr22.
