# Phase8G pruned 5-feature candidate v1

Historical note: this precision candidate is not the current repository default.
The current default is tracked in `docs/current_default.md`.

## Status

Interpretable precision candidate.

Balanced/default branch remains:
- unpruned Phase8F threshold 0.50

Recommended pruned modes:
- Phase8G-pruned threshold 0.80 = precision mode
- Phase8G-pruned threshold 0.90 = high-precision mode

## Selected features

- panel_margin
- support_conflict_ratio
- best_vs_second_margin
- phase8e_unweighted_plus_support
- phase8e_unweighted_confidence

## Interpretation

The pruned model reduces the feature set from 33 usable features to 5 while slightly improving candidate-level CV AUC/AP.

At downstream phasing level, the pruned model is most useful as a stricter precision gate. At threshold 0.80 it reduces hamming and switch errors relative to unpruned Phase8F t0.80, with a small truth_correct_pct decrease.
