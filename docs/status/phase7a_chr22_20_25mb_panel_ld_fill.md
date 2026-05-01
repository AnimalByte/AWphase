# Phase7A chr22 20-25Mb panel-LD fill status

## Best balanced candidate

Candidate: d20000_s20_c085_m5_b2

Compared against Phase6C 35x and real WhatsHap Illumina 35x on HG002 chr22:20-25Mb.

## Metrics

Phase7A balanced:
- n_pred_sites_nonzero: 3391
- n_exact_overlap_sites_phased: 2398
- hamming_errors: 7
- hamming_error_rate: 0.0029190992493744786
- switch_errors: 8
- switch_error_rate: 0.004308023694130318
- phased_site_accuracy_pct: 99.70809007506254
- truth_correct_pct: 67.50423489553924
- raw_block_n50_bp: 1336694
- max_block_span_bp: 1336694
- median_block_span_bp: 764

WhatsHap 35x:
- n_pred_sites_nonzero: 3370
- n_exact_overlap_sites_phased: 2369
- hamming_errors: 23
- hamming_error_rate: 0.009708737864077669
- switch_errors: 17
- switch_error_rate: 0.009675583380762664
- phased_site_accuracy_pct: 99.02912621359224
- truth_correct_pct: 66.23376623376623
- raw_block_n50_bp: 1336694

## Interpretation

Phase7A balanced improves recall beyond WhatsHap on this benchmark window while retaining lower hamming/switch error and the same long-block N50.

Phase6C remains the trusted read-backed wMEC core. Phase7A is a panel-LD fill layer for sites that read evidence did not phase directly.

## Block-mixing vs WhatsHap

- blocks_evaluated: 678
- total_overlap_sites: 3245
- residual_error_rate_after_best_block_orientation: 0.02403697996918336
- adjacent_pair_mismatch_rate: 0.020646669263731984
- mixed_blocks: 45

This is a disagreement-against-WhatsHap diagnostic, not direct truth error. It should be used to target calibration and panel-fill risk control.
