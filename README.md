# AWPhase

AWPhase is a research-stage, panel-assisted phasing project for HG002 short-read
phasing experiments. The current project target is Phase9D.

## Quick Start

Install the managed environment:

```bash
pixi install
```

Stage the reusable runtime data needed by the phasing scripts: HGDP + 1KGP
SHAPEIT5 reference panels, Beagle GRCh38 genetic maps, panel metadata, and the
GRCh38 no-alt reference FASTA:

```bash
CHROMS="chr1 chr20 chr22 chr6" THREADS=12 bash scripts/setup_awphase_runtime_data.sh
```

This does not run a benchmark and does not download HG002 read BAMs or GIAB
truth data. To stage only one chromosome:

```bash
CHROMS="chr6" THREADS=12 bash scripts/setup_awphase_runtime_data.sh
```

Raw genomics data is intentionally kept out of git. Runtime data is written under
`data/` and `references/`, both ignored by `.gitignore`.

## Current Default

The current documented default is:

```text
Phase9D default:
  base method: Phase8F_ancestry_weighted_multiwindow
  threshold:   0.80
```

See [docs/current_default.md](docs/current_default.md) for the promoted modes,
important caveats, and the distinction between the promoted Phase9D default and
the experimental Phase9 bridge-forest scripts.

## Repository Layout

- `rust/`: Rust phasing core and CLI prototype.
- `python/awphase_py/`: Python experiment, training, feature, and evaluation tools.
- `scripts/`: shell and Python runners for benchmark windows and phase experiments.
- `configs/`: benchmark configs and current profile metadata.
- `docs/`: architecture, data, benchmark, and status notes.
- `freezes/`: frozen models, metrics, and summary artifacts.
- `results/`: retained summary outputs and phase comparison tables.

## Export Caveat

This checkout was restored from a portable export. Raw data and large external
artifacts are intentionally not included: `data/`, `baselines/`, `references/`,
raw BAM/CRAM/VCF/BCF/FASTQ/reference files, and large intermediates must be
staged separately before full benchmarks can run.

## Environment

Use [pixi.toml](pixi.toml) as the primary environment definition. The older
[environment.yml](environment.yml) is retained as a conda-compatible reference.
The Python tools are currently run directly with `PYTHONPATH=python`; they are
not yet packaged as an installed command.

## Lightweight Checks

Install the environment and validate the current Phase9D profile with:

```bash
pixi install
pixi run validate-phase9d
pixi run validate-benchmark-bundles
```

Inspect a retained Phase9D benchmark window without downloading data:

```bash
pixi run phase9d-window list --split train
pixi run phase9d-window plan chr22_20_25mb
pixi run phase9d-window check chr22_20_25mb
```

`plan` prints the manifest-derived commands and expected outputs. `check`
requires the external BAM, BCF, VCF, BED, and reference FASTA paths to already
exist locally.

Run the default lightweight project checks with:

```bash
pixi run check
```

Run the full local CI set with:

```bash
pixi run ci
```

## HG002 Benchmark Setup

For the original chr1/chr20/chr22 HG002 benchmark bundle, use:

```bash
THREADS=12 bash scripts/setup_awphase_benchmark.sh setup
```

For the chr6 MHC hard-region benchmark used in the latest PBWT/HMM comparison,
use:

```bash
THREADS=12 bash scripts/setup_chr6_mhc_benchmark.sh
```

That chr6 script downloads `chr6:25,000,000-34,999,999` slices first and then
downsamples the local slices to 30x for both Illumina 151bp and PacBio Revio
HiFi. It stages data only; benchmarks are run separately from the manifest.

The retained chr6 summary tables are:

- `results/phase8/pbwt_hmm_chr6_mhc_manifest/phase8_pbwt_hmm.hard_chr6.average.tsv`
- `results/phase8/pbwt_hmm_chr6_mhc_manifest/phase8_pbwt_hmm.hard_chr6.per_window.tsv`

Panel-reliance diagnostics are retained separately:

- `results/phase8/panel_reliance_chr6_mhc_manifest/panel_reliance.hard_chr6.average.tsv`
- `results/phase8/panel_reliance_manifest/panel_reliance.train.average.tsv`

The Phase8 PBWT/HMM runner now emits the original PBWT/HMM paths, v2
exact-prefix PBWT with forward-backward HMM, and an experimental v3
bidirectional-prefix PBWT path. On the retained chr6 MHC run, v3 matched v2 but
was slower, so v2 remains the better current default.

To run AWPhase itself on PacBio Revio HiFi 30x input instead of the manifest's
Illumina BAM, use:

```bash
AWPHASE_READSET=pacbio_hifi30x \
  bash scripts/run_phase8_pbwt_hmm_manifest.sh \
  results/phase8f_manifests/chr6_mhc_windows.window_beds.tsv hard_chr6
```

The default PacBio runner profile uses smaller Phase6C WMEC chunks
(`AWPHASE_MAX_COMPONENT_SITES=64`, `AWPHASE_LOCAL_REFINE_ITERS=20`) and writes
to `results/phase8/pbwt_hmm_pacbio_hifi_mcs64_manifest/`,
`results/phase8/panel_reliance_pacbio_hifi_mcs64_manifest/`, and
`results/phase8/confidence_calibration_pacbio_hifi_mcs64_manifest/`. On chr6
MHC this reduced PacBio-AWPhase PBWT+HMM v2 hamming error from about 20.1% in
the original large-component diagnostic profile to about 1.86%, with 93.29%
phased coverage and 91.55% truth-correct coverage. It also improves on the
mcs128 profile, which was about 4.48% hamming error at 93.35% phased coverage.
It is still below WhatsHap PacBio HiFi, so PacBio input mode remains diagnostic
rather than promoted.

Confidence calibration summaries are produced by
`scripts/phase8/calibrate_phase8_confidence_v1.py`. The retained chr6 MHC
PacBio mcs64 calibration shows the remaining AWPhase PacBio PBWT+HMM v2 calls
are overconfident by about 1.1 percentage points, while the mcs128 profile was
overconfident by about 3-4 percentage points and the original PacBio
large-component profile was overconfident by about 18-19 percentage points.
