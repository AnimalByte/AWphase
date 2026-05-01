# AWPhase Phase8 implementation plan

Historical note: this Phase8 plan records the stable branch at the time it was
written. The current repository default is tracked in `docs/current_default.md`.

## Historical stable branch

Stable branch at the time:

- Phase7G donor XGBoost
- Default threshold: 0.50
- Precision threshold: 0.80

Phase7G is a high-precision, moderate-recall, short-read panel-assisted phaser.

Six-window benchmark:

- chr1_15_20mb
- chr20_25_30mb
- chr20_30_35mb
- chr20_45_50mb
- chr22_15_20mb
- chr22_20_25mb

Six-window average at the time:

AWPhase Phase7G default:
- weighted_truth_correct_pct: 70.23
- weighted_hamming_error_rate: 0.00477
- weighted_switch_error_rate: 0.00255

WhatsHap Illumina 30x:
- weighted_truth_correct_pct: 83.51
- weighted_hamming_error_rate: 0.02649
- weighted_switch_error_rate: 0.02032

WhatsHap GIAB CCS/HiFi 30x:
- weighted_truth_correct_pct: 97.48
- weighted_hamming_error_rate: 0.02208
- weighted_switch_error_rate: 0.00719

## Phase8 goals

Phase8 moves from PBWT/HMM-inspired features to actual PBWT/HMM-style panel donor scoring.

Excluded from Phase8:
- trio/pedigree phasing logic
- Nextflow production pipeline

External holdout and full chromosome benchmarking are intentionally delayed until the engine is more complete.

## Phase8A: true PBWT donor retrieval

Goal:
- Build a real positional PBWT-style donor retrieval layer from HGDP+1KGP phased panel haplotypes.
- Given Phase6C anchor states around a candidate site, retrieve locally matching panel haplotypes efficiently.

Deliverables:
- PBWT panel index builder
- per-window PBWT donor retrieval features
- donor neighbor counts
- donor top-k match features
- serialized index/cache per chromosome/window

Success criteria:
- PBWT retrieval features are non-empty across all six windows.
- PBWT features improve candidate-level CV or final gated metrics versus Phase7G.

## Phase8B: genetic-map HMM donor scoring

Goal:
- Replace bp-distance proxy with genetic-map-aware recombination distances.
- Score donor haplotype paths using an HMM-like model.

Deliverables:
- genetic map loader
- bp to cM interpolation
- HMM transition probability from genetic distance
- emission model from anchor agreement/disagreement
- HMM log-likelihood margin features

Success criteria:
- HMM features rank among important features.
- Error rate improves at similar truth_correct_pct.

## Phase8C: Beagle-style haplotype copying

Goal:
- Add haplotype copying features inspired by Beagle.
- Use top local donor haplotypes as a compact state space.
- Estimate candidate phase from donor copying probabilities.

Deliverables:
- top-k donor state compression
- donor posterior entropy
- effective donor count
- left/right donor consistency
- donor copying posterior features

Success criteria:
- Improves recall without increasing hamming/switch materially.

## Phase8D: SHAPEIT5-style rare-variant projection

Goal:
- Treat Phase6C common/strong anchors as a scaffold.
- Project low-support/rare candidate sites onto the scaffold using panel donor evidence and read-backed constraints.

Deliverables:
- common/scaffold site classifier
- rare/fill candidate classifier
- projection confidence features
- scaffold-aware candidate grouping

Success criteria:
- Improves fill recall while preserving Phase7G-level error.

## Phase8E: panel ancestry weighting / ELAI-inspired priors

Goal:
- Weight donor haplotypes by ancestry compatibility.
- Start with panel metadata/global ancestry.
- Later add local ancestry priors if ELAI/RFMix/gnomix outputs exist.

Deliverables:
- panel sample ancestry metadata loader
- global ancestry prior weights
- local ancestry prior interface
- ancestry-weighted donor support/margin features

Success criteria:
- Weighted donor features improve over unweighted donor features.
- No loss of precision in admixed-like windows.

## Phase8F: larger multi-window learned model

Goal:
- Train a larger Phase8 model across more windows and chromosomes.
- Avoid leakage by grouping on physical windows.

Deliverables:
- selected training windows manifest
- automated run matrix
- combined training table
- grouped cross-validation
- threshold selection report

Success criteria:
- Stable grouped CV.
- Default threshold selected by final phasing metrics, not candidate AUC only.

## Phase8G: clean AWPhase CLI

Goal:
- Provide one command entry point for window-level AWPhase.

Example target CLI:

awphase phase-window \
  --chrom chr20 \
  --start 30000000 \
  --end 34999999 \
  --bam sample.bam \
  --panel-bcf panel.bcf \
  --variant-json variants.json \
  --truth-vcf truth.vcf.gz \
  --truth-bed truth.bed \
  --mode phase8

Deliverables:
- awphase CLI module
- subcommands for build-fragments, solve, panel-fill, donor-features, train, evaluate
- consistent output directories

Success criteria:
- One command can reproduce a window run.

## Phase8H: manuscript-ready tables/plots

Goal:
- Generate reproducible benchmark tables and plots.

Deliverables:
- per-window comparison tables
- average comparison tables
- precision/recall/error tradeoff plots
- truth_correct versus error plots
- block N50 plots
- feature importance plots
- model summary table

Success criteria:
- Figures can be dropped into a report/paper.

## Phase8I: external holdout benchmark

Goal:
- Test Phase8 on unseen windows not used in training.

Candidate holdout windows:
- chr1_100_105mb
- chr1_165_170mb
- chr22_25_30mb
- chr22_30_35mb
- chr22_35_40mb

Success criteria:
- Phase8 default generalizes beyond the original six windows.

## Phase8J: full chromosome benchmark

Goal:
- Run chromosome-scale benchmark after holdout validation.

Candidate chromosomes:
- chr20 first
- chr22 second
- chr1 last

Success criteria:
- Scales without pathological memory/runtime.
- Maintains precision advantage.
- Produces credible full-chromosome comparison against WhatsHap Illumina and GIAB CCS/HiFi.
