# Phase7E v4b curated 30x stable

Stable/default learned gate for Phase7A panel-fill candidates.

Training:
- 6 physical 30x windows
- chr1, chr20, chr22
- 2,451 candidate rows
- 2,149 positives
- 302 negatives

Default:
- Phase7E-balanced threshold = 0.70

Alternative:
- Phase7E-precision threshold = 0.80

Interpretation:
v4b is the current stable default because it gives a consistent precision/recall tradeoff and is less feature-heavy than Phase7F.
