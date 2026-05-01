# Phase8A/B PBWT-HMM experimental v1

Historical note: this document describes an experimental branch and is not the
current repository default. The current default is tracked in
`docs/current_default.md`.

## Status

Experimental. Do not promote over Phase7G.

## What was implemented

This branch added PBWT/HMM-style donor scoring features on top of the Phase7G candidate table:

- genetic-map-aware anchor selection
- cM-distance decay
- top donor haplotype support
- PBWT-like plus/minus support
- HMM-like log margin
- PBWT confidence/entropy/effective donor features

## Important caveat

This is not yet a true Durbin PBWT index implementation. It is a PBWT/HMM-inspired donor scoring layer using panel haplotypes and genetic map distance.

## Result

Phase8A/B v1 was close to Phase7G but did not clearly outperform it.

Stable branch at the time remained:

- Phase7G donor XGBoost
- default threshold: 0.50
- precision threshold: 0.80

## Interpretation

The Phase8A/B features are useful, but this v1 implementation is not enough to replace Phase7G. Next work should focus on:

1. actual PBWT prefix/divergence indexing, or
2. SHAPEIT5-style scaffold/rare projection, or
3. ancestry-weighted donor scoring.
