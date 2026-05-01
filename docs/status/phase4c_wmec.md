# Phase 4c wMEC patch

Fixes Phase 4b evidence starvation by:
- building regions through the full read-connectivity graph, not candidate-only edges
- relaxing but still filtering allele observations
- merging reads only under stronger same-haplotype log-odds
- increasing coverage cap and read cap
- adding per-region diagnostics
- making the wMEC solve dominant while treating ranker/donor/scaffold as penalties
