# AWPhase portable export

This folder contains the AWPhase code, scripts, configs, documentation, frozen models, and summary metrics.

Excluded intentionally:
- data/
- baselines/
- references/
- raw BAM/CRAM/VCF/BCF/FASTQ/reference files
- large intermediate files such as fragments, local calls, and full window outputs

Important folders:
- python/awphase_py: Python implementation
- scripts/: run scripts and phase scripts
- docs/status/: status notes and phase decisions
- freezes/: frozen candidate models and summaries
- results/phase8: Phase8F/Phase8G metrics and pruning summaries
- results/phase9: Phase9 bridge/forest summaries
- results/phase8f_manifests: benchmark window manifest
