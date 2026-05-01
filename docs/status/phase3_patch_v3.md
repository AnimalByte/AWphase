# Phase 3 patch v3

This patch applies the Phase 3 upgrades in one shot:
- local group-conditioned donor prior aligned to `hg002_chr20_target_sites.tsv`
- explicit reference-context feature extraction from FASTA
- scaffold-aware candidate features with anchor agreement and switch-risk fields
- XGBoost ranking patch using query-group ranking
- continuity-aware candidate application with margin/anchor/donor guards

Expected behavior:
- donor prior summary becomes non-zero
- candidate features gain anchor and donor fields
- ranker output becomes more conservative than the previous free-flip version
- switch should be protected relative to the prior Phase 3 run
