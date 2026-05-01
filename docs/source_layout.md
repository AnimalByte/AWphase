# Source Layout And Artifact Policy

This repository was restored from a portable export. Treat the extracted source
tree as the active source of truth and the export tarball as an archival
handoff artifact, not as the development unit.

## Track Normally

- Rust crates under `rust/` plus root `Cargo.toml` and `Cargo.lock`
- Python source under `python/awphase_py/`
- small runner and validator scripts under `scripts/`
- Pixi and conda environment definitions
- configs, docs, tests, workflow stubs, and curated small summary tables

## Do Not Track By Default

- external data roots: `data/`, `baselines/`, `references/`
- BAM, CRAM, FASTQ, VCF/BCF, FASTA, and their indexes
- generated run directories such as `results/runs/`, `results/experiments/`,
  `results/comparisons/`, and `results/phase7a_windows/`
- Python caches, Rust `target/`, Pixi environment directories, logs, and scratch
  output

## Curated Results

Small result summaries in `results/` and `freezes/` can be tracked when they are
part of a documented benchmark freeze or decision note. Large intermediates
should stay local and be recreated from the manifest plus staged external data.

## Current Entry Points

- `pixi run ci`: local source validation
- `pixi run phase9d-window plan <label>`: manifest-derived Phase9D window plan
- `pixi run phase9d-window check <label>`: local data presence and bundle check

None of these commands downloads genomics data.
