# AWPhase Data Contract

The portable export contains source code, scripts, configs, documentation,
frozen models, and retained summary metrics. It intentionally does not contain
large genomics inputs or generated intermediates.

## Excluded External Roots

The following roots are expected to be staged separately:
- `data/`
- `baselines/`
- `references/`

Common excluded file types include BAM/CRAM, VCF/BCF, FASTQ, FASTA, index files,
fragment matrices, full local-call outputs, and large benchmark intermediates.

## Required Data Classes

A full Phase9D benchmark requires:
- target sample read alignments, usually HG002 Illumina 30x windows
- target variant JSON derived from the target VCF
- phased panel BCF/VCF files, such as HGDP+1KGP panel slices
- truth VCF and matching benchmark BED files
- reference FASTA and indexes when extracting read observations
- retained Phase8F/Phase9 feature and prediction tables for bridge experiments

## Bundle Rules

Truth VCF and benchmark BED files must come from the same benchmark release
bundle. Do not mix v5.0q truth VCFs with older v4.2.1 BEDs or vice versa.

Panel data must exclude exact target-sample leakage. For HG002, the leakage check
should cover common aliases such as HG002, NA24385, and GM24385.

## Manifest Direction

The next reproducibility step is a checked data manifest with:
- local path
- source URL or accession
- checksum
- required indexes
- benchmark release version
- sample exclusions
- license or access notes
