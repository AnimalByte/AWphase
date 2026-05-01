# AWPhase Architecture

AWPhase is split into a stable engine surface and a research workflow surface.

## Stable Engine Direction

The stable path should live in Rust:
- configuration and manifest validation
- deterministic panel and read-backed phasing logic
- local decoding, block construction, bridge scoring, and stitching
- stable JSON/TSV output schemas
- benchmark and data-contract validators

The current Rust workspace contains:
- `rust/awphase-core`: core types, candidate scoring, local decode, bridge, IO,
  reports, and summaries
- `rust/awphase-cli`: prototype CLI with `candidates`, `local`, `bridge`, and
  `benchmark` subcommands

## Research And ML Surface

Python remains the right place for fast-moving work:
- XGBoost training and scoring experiments
- pandas/sklearn feature tables and grouped cross-validation
- comparison summaries and manuscript tables
- BAM/VCF extraction through `pysam`
- exploratory phase-specific scripts

The Python tools currently run through `PYTHONPATH=python python ...` rather
than an installed package.

## Current Target

The current project target is Phase9D. The promoted default is documented in
`docs/current_default.md` and summarized by
`configs/current_phase9d.profile.yaml`.

## Known Gaps

- The Nextflow placeholder at `workflow/main.nf` is empty.
- The portable export excludes raw data and large intermediates.
- The stable one-command Phase9D workflow has not yet been promoted into a
  single CLI command.
- Several historical status documents describe older defaults. Treat
  `docs/current_default.md` as the current pointer.
