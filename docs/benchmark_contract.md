# AWPhase Benchmark Contract

This contract defines how AWPhase comparisons should be made before promoting a
new default.

## Current Baseline

The current promoted AWPhase mode is Phase9D default:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.80
```

Experimental Phase9 bridge-forest outputs may be reported, but they should not
replace the promoted default unless a later decision document says so.

## Required Comparisons

Promoted modes should be compared against:
- WhatsHap Illumina baseline
- PacBio/HiFi-backed WhatsHap benchmark where available
- unbridged Phase8F thresholds `0.70`, `0.80`, and `0.90`
- the previous promoted AWPhase default

## Required Metrics

Report at minimum:
- hamming errors
- hamming error rate
- switch errors
- switch error rate
- phased-site accuracy
- truth-correct percentage
- block N50
- median block span
- maximum block span
- number of bridged windows or accepted bridge edges where applicable
- windows improved versus degraded

## Validity Rules

- Use grouped window splits for learned models.
- Keep held-out windows separate from training windows.
- Do not mix truth VCF and BED benchmark bundle versions.
- Record the exact threshold and model artifact used.
- Record missing windows rather than silently dropping them from averages.
- Keep comparison sets explicit: all available, common with Illumina, common with
  HiFi/PacBio, and common across all methods are different claims.

Validate static benchmark bundle references with:

```bash
pixi run validate-benchmark-bundles
```
