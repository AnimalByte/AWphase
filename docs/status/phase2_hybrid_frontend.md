# Phase 2 hybrid-read-frontend closeout

## Promoted branches
- default local-decode branch: `selected`
- confidence / arbitration branch: `tagged`
- continuity / bridge branch: `superreads`

## Experimental branch
- `hybrid` remains experimental and is not promoted

## Benchmark policy
- HG002 primary benchmark: `v5.0q` where matching local truth VCF + benchmark BED are staged
- HG002 secondary benchmark: `v4.2.1` for comparability and for chromosomes/intervals where local v5.0q benchmark assets are not yet staged

## Outputs
- `freezes/phase2_hybrid_frontend/summary.json`
- `results/metrics/phase2_hybrid_frontend_summary.tsv`
- `results/metrics/v5q_chr20/*.metrics.json`
- `results/metrics/v5q_asset_audit.tsv`
