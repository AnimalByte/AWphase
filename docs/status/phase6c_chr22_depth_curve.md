# Phase6C chr22 20-25Mb depth curve

## Summary

Phase6C was compared against WhatsHap on the same downsampled chr22 20-25Mb BAMs from HG002 Illumina 35x.

## Main result

Across 20x, 25x, 30x, and 35x, Phase6C has lower hamming/switch error than WhatsHap at the same depth, but WhatsHap phases more sites and has higher truth_correct_pct.

## Interpretation

Phase6C is currently a high-precision, lower-recall phaser. WhatsHap is higher-recall and slightly lower precision on this benchmark.

## Depth curve

20x:
- Phase6C hamming: 0.0
- Phase6C switch: 0.0
- Phase6C truth_correct_pct: 59.93788819875776
- Phase6C raw_block_n50_bp: 1336002
- WhatsHap hamming: 0.005620406398616515
- WhatsHap switch: 0.006630500301386378
- WhatsHap truth_correct_pct: 64.93506493506493
- WhatsHap raw_block_n50_bp: 1335391

25x:
- Phase6C hamming: 0.0004616805170821791
- Phase6C switch: 0.0006203473945409429
- Phase6C truth_correct_pct: 61.12365894974591
- Phase6C raw_block_n50_bp: 1336002
- WhatsHap hamming: 0.005972696245733789
- WhatsHap switch: 0.006462984723854289
- WhatsHap truth_correct_pct: 65.78204404291361
- WhatsHap raw_block_n50_bp: 1336002

30x:
- Phase6C hamming: 0.00045641259698767686
- Phase6C switch: 0.0006079027355623101
- Phase6C truth_correct_pct: 61.82947487295314
- Phase6C raw_block_n50_bp: 1336002
- WhatsHap hamming: 0.007627118644067797
- WhatsHap switch: 0.008083140877598153
- WhatsHap truth_correct_pct: 66.12083568605308
- WhatsHap raw_block_n50_bp: 1336002

35x:
- Phase6C hamming: 0.00045310376076121433
- Phase6C switch: 0.0005995203836930455
- Phase6C truth_correct_pct: 62.28119706380576
- Phase6C raw_block_n50_bp: 1336694
- WhatsHap hamming: 0.009708737864077669
- WhatsHap switch: 0.009675583380762664
- WhatsHap truth_correct_pct: 66.23376623376623
- WhatsHap raw_block_n50_bp: 1336694

## Caveat

The earlier original 30x BAM gave Phase6C a much smaller N50. The downsampled 30x from the 35x BAM gives long blocks. This suggests the original 30x BAM differs in effective interval coverage/connectivity and should be diagnosed separately.

## Recall-gap diagnostic

The Phase6C-vs-WhatsHap recall gap on chr22 20-25Mb is not primarily long-block recall.

WhatsHap-only sites relative to Phase6C 35x:
- total: 218
- SNP: 218
- singleton WhatsHap blocks: 138
- WhatsHap blocks >100 sites: 0

Interpretation:
WhatsHap's additional sites are mostly singleton/tiny-block SNP calls. Phase6C is not missing many sites inside the long WhatsHap block. The remaining gap in truth_correct_pct is therefore partly driven by tiny-block/singleton emissions rather than long-range haplotype phasing.

Phase6D same-template SNP fill produced zero accepted fills, suggesting that SNP sites directly attachable to the Phase6C backbone were already captured by the Phase6C BAM-template wMEC graph.

## Panel leakage check

The HGDP+1KGP chr22 panel sample list contains sample IDs such as HG00233, HG00234, etc. These are not the GIAB HG002 sample. Exact checks should be used instead of broad prefix matching.

Required exact leakage check:
- HG002
- NA24385
- GM24385

If exact matching returns no samples, the panel is acceptable for Phase7A testing.
