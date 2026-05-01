# Phase 5b chr22 20-25Mb status

## Valid comparison

Real chr22 Illumina WhatsHap baselines were generated using chr22 interval BAMs. Earlier zero-call WhatsHap rows were invalid because chr20 BAMs were accidentally used.

## Current chr22 AWPhase result

AWPhase readbridge-superreads improves over AWPhase tagged baseline:

- tagged hamming: 0.20992761116856257
- readbridge-superreads d1000 hamming: 0.20367816091954022
- tagged switch: 0.272608125819135
- readbridge-superreads d2500 switch: 0.2475975127190503

## External WhatsHap gap

WhatsHap Illumina 30x is dramatically better:

- hamming: 0.008959044368600682
- switch: 0.009313154831199068
- phased-site accuracy: 99.10409556313992

## Diagnostic conclusion

Site-matched and block-mixing diagnostics show AWPhase is not mainly losing due to fewer sites or whole-block orientation. The major failure is mixed/noisy within-block phase patterns.

Block-mixing summary against WhatsHap:

- blocks evaluated: 435
- overlap sites: 2672
- residual errors after best block orientation: 630
- residual error rate after best block orientation: 0.23577844311377247
- adjacent-pair mismatch rate: 0.25435851586946806
- mixed blocks: 241

## Phase 6 direction

Replace the graft-centered decision layer with a true read-fragment weighted MEC core:

1. Build strict allele observation matrix.
2. Construct read-fragment connected components.
3. Solve weighted MEC per component.
4. Emit coherent PS/block IDs directly from components.
5. Add scaffold/donor priors only as soft terms after the read-backed MEC core exists.
