# Phase7G donor XGBoost stable

Phase7G is the current stable/default learned gate for AWPhase panel-assisted short-read phasing.

## Concept

Phase7G combines:

- Phase6C read-backed wMEC backbone
- Phase7A panel-fill candidates
- Phase7F anchor/block geometry
- PBWT/HMM-inspired donor-copying features
- XGBoost accept/reject gate

## Training/evaluation windows

Six 30x windows:

- chr1_15_20mb
- chr20_25_30mb
- chr20_30_35mb
- chr20_45_50mb
- chr22_15_20mb
- chr22_20_25mb

## Recommended modes

Default:
- Phase7G donor gate
- threshold = 0.50

Precision mode:
- Phase7G donor gate
- threshold = 0.80

## Six-window average comparison

AWPhase Phase7G default t0.50:

- weighted_truth_correct_pct: 70.23
- weighted_hamming_error_rate: 0.00477
- weighted_switch_error_rate: 0.00255

WhatsHap Illumina 30x:

- weighted_truth_correct_pct: 83.51
- weighted_hamming_error_rate: 0.02649
- weighted_switch_error_rate: 0.02032

WhatsHap GIAB CCS/HiFi 30x:

- weighted_truth_correct_pct: 97.48
- weighted_hamming_error_rate: 0.02208
- weighted_switch_error_rate: 0.00719
- mean_raw_block_n50_bp: 556,146 bp

## Interpretation

AWPhase Phase7G is a high-precision, moderate-recall, short-read panel-assisted phaser.

Compared with WhatsHap Illumina 30x, AWPhase Phase7G has much lower hamming and switch error, but lower truth-site recall.

Compared with GIAB CCS/HiFi WhatsHap 30x, AWPhase Phase7G remains much lower recall and shorter-block, while HiFi WhatsHap represents the long-read upper-bound baseline.
