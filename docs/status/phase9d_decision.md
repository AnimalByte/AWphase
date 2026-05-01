# AWPhase Phase9D Decision

Current repository pointer: `docs/current_default.md`.

## Summary

Phase9D promotes the Phase8F ancestry-weighted multiwindow method as the current AWPhase default path.

The best practical default is:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.80
```

This mode gives the best balance between improved phase continuity and controlled error risk.

## Promoted AWPhase modes

### 1. Default balanced mode

Use:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.80
```

This is the recommended default AWPhase mode for Phase9D.

### 2. Conservative mode

Use:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.90
```

Use this when accuracy is more important than block extension, especially for rare-disease interpretation where false bridges could affect cis/trans or compound-heterozygous interpretation.

### 3. Exploratory bridge-extension mode

Use:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.70
```

Use this only for exploratory analysis where larger blocks are useful, but where a higher risk of incorrect joins is acceptable.

## Decision

Phase9D should move forward with:

```text
Default:       threshold 0.80
Conservative:  threshold 0.90
Exploratory:   threshold 0.70
```

The production-facing AWPhase default should be threshold `0.80`.

The conservative rare-disease interpretation mode should be threshold `0.90`.

The exploratory bridge-extension mode should be threshold `0.70`.

## Interpretation relative to WhatsHap and PacBio

AWPhase should not be treated as replacing PacBio or true long-read phasing.

Instead:

- PacBio / HiFi remains the upper benchmark for long-range phase truth.
- WhatsHap remains the read-backed short-read/local phasing baseline.
- AWPhase is a statistical ancestry-aware bridge-extension layer on top of short-read data.
- The goal of AWPhase is to extend useful phase information beyond local read-backed blocks while avoiding unsafe over-bridging.

A successful AWPhase result should improve short-read phasing utility without moving too far away from the PacBio/HiFi truth benchmark.

## Next benchmark comparison

The promoted Phase9D modes should be compared against:

1. WhatsHap Illumina baseline
2. PacBio/HiFi-backed phasing benchmark
3. Original unbridged Phase8F baseline
4. Phase8F thresholds 0.70, 0.80, and 0.90

Key metrics to report:

- hamming errors
- switch errors
- phased-site accuracy
- truth-correct percentage
- block N50
- median block span
- number of bridged windows
- number of windows where bridging improved vs degraded

## Final Phase9D call

Promote:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.80
```

as the AWPhase Phase9D default.

Keep `threshold 0.90` as the conservative interpretation mode.

Keep `threshold 0.70` as exploratory only.
