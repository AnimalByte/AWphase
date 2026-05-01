# Phase 3 internal scaffold and candidates patch v2

This patch improves the first internal Phase 3 pass by:
- increasing scaffold coverage via obs-level support integration
- making donor prior extraction schema-tolerant and position-aware
- adding explicit reference-context features from FASTA + variant JSON
- auto-detecting the phased-truth sample during label generation
- wiring the reference-context stage into the Phase 3 chr20 internal runner
