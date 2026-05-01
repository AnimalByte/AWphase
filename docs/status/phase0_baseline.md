# Phase 0 baseline freeze

## Decisions
- default local decode branch: `selected`
- confidence / arbitration branch: `tagged`
- continuity / bridge branch: `superreads`

## Why
These branches generalized across chr1, chr20, and chr22 better than the plain V2 read path.

## Files
- `freezes/phase0_baseline/summary.json`
- `results/metrics/phase0_branch_summary.tsv`
