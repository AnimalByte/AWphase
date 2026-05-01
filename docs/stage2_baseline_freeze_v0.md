# awphase Stage-2 Baseline Freeze v0

## What is frozen
This freeze captures the first real Stage-2 local structured decoder baseline.

### Decoder characteristics
- two-pass local decode flow
- candidate-conditioned local refinement
- explicit window-level candidate states
- window-level Viterbi-style donor-path decoding
- posterior-like window margins
- confidence-aware abstention
- path-vs-posterior agreement reporting

### What is NOT yet included
- PBWT-backed candidate retrieval
- composite / compressed donor states as real core states
- full Li-Stephens / HMM copy model
- mature bridge layer
- truth-calibrated thresholds

## Why frozen now
This is the first version that:
- has a coherent best donor-copy path per window
- preserves abstention when posterior support is weak
- is stable enough to benchmark on HG002 chr20 before PBWT changes it

## Next required task
Run HG002 chr20 window-level truth benchmarking on this frozen baseline before changing retrieval or decoder structure.
