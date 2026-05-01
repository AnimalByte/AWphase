# AWPhase Current Default

## Canonical Target

The current AWPhase target is Phase9D.

The promoted Phase9D default is:

```text
Phase8F_ancestry_weighted_multiwindow threshold 0.80
```

This is the production-facing default from the Phase9D decision note. It should
be treated as the baseline current mode until a newer decision document replaces
it.

## Modes

Default:
- method: `Phase8F_ancestry_weighted_multiwindow`
- threshold: `0.80`
- use when balancing useful phase continuity against controlled error risk

Conservative:
- method: `Phase8F_ancestry_weighted_multiwindow`
- threshold: `0.90`
- use when false joins are more costly, such as rare-disease interpretation

Exploratory:
- method: `Phase8F_ancestry_weighted_multiwindow`
- threshold: `0.70`
- use only when larger blocks are useful and higher join risk is acceptable

## Phase9 Bridge-Forest Status

The Phase9D bridge-forest/XGBoost utilities exist in the codebase, but they are
not the promoted default in the current decision note. Treat them as experimental
bridge-extension work unless a later benchmark decision promotes them.

Relevant files:
- `python/awphase_py/phase9/apply_phase9d_utility_bridge_forest_v1.py`
- `python/awphase_py/phase9/train_phase9d_wide_bridge_xgb_noleak_safe.py`
- `results/phase9/phase9d_wide_bridge_xgb.noleak.metrics.json`
- `results/phase9/phase9d_utility_forest_on_phase8f_t080.average.tsv`

## Window Planning

Use the manifest-driven wrapper to inspect or validate a benchmark window without
downloading data:

```bash
pixi run phase9d-window plan chr22_20_25mb
pixi run phase9d-window check chr22_20_25mb
```

`plan` is safe before data staging because it only reads repository metadata.
`check` fails until the external BAM, panel, truth, and reference files are
present locally.

## Benchmark Caveat

AWPhase is a statistical ancestry-aware bridge-extension layer on top of
short-read phasing. It should not be described as replacing PacBio, HiFi, or
other true long-read phasing. WhatsHap remains the read-backed local phasing
baseline, and PacBio/HiFi remains the long-range truth-oriented upper benchmark.

## Source Of Truth

This summary is a stable pointer for the repository. The detailed decision lives
in `docs/status/phase9d_decision.md`.

The machine-readable profile is `configs/current_phase9d.profile.yaml`. With the
Pixi environment installed, validate it with:

```bash
pixi run validate-phase9d
```
