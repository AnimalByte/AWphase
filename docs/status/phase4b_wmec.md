# Phase 4b: WhatsHap-inspired bounded-coverage wMEC layer

This patch adds six WhatsHap-inspired mechanics to AWPhase:

1. Region decision core is a wMEC-style weighted block objective
2. Active coverage is bounded before solving
3. Regions are built from read-connectivity components, then split by scaffold/size
4. The matrix excludes weak/ambiguous observations
5. Scaffold blocks are converted into synthetic reads inside the same optimization
6. Read merging is probabilistic, based on same-vs-different haplotype likelihood

Implementation notes:
- Exact enumeration is used on small regions (bounded by `max_exact_sites`)
- Larger regions use a beam-search fallback
- Donor/ranker/anchor terms are used as priors/penalties, not as the primary objective
