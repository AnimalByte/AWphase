# Phase 4 MEC-like resolver

This patch adds the next missing pieces modeled after WhatsHap/SHAPEIT5/Beagle ideas while keeping AWPhase internal:

- scaffold blocks as synthetic reads
- local donor prior aligned to target-site ordering
- haplotag-aware realign-like support features
- scaffold-aware candidate features v2
- XGBoost rank:pairwise candidate ranker
- region-level MEC-like resolver over small ambiguous windows
- conservative application with margin + anchor/donor guards

This sits on top of the existing Phase 3 internal pipeline.
