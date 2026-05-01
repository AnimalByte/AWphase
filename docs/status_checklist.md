# awphase status checklist

## Frozen-plan alignment

### Stage 0: architecture freeze
- [x] Frozen math backbone
- [x] Frozen execution plan
- [x] Module boundaries exist

### Stage 1: candidate engine
- [x] Candidate-state interface
- [x] Panel haplotype loader
- [x] Deterministic bounded candidate ranking
- [x] Ancestry-weighted donor scoring
- [x] Simple composite state output
- [ ] Real PBWT-backed retrieval
- [ ] Real donor compression / composite logic

### Stage 2: local decoder
- [x] Two-pass decode flow
- [x] Candidate-conditioned local refinement
- [x] Explicit local window states
- [x] Confidence / abstention output
- [ ] Reduced-state Li-Stephens / HMM local decoder
- [ ] Explicit small window-state decoder with competing retained states
- [ ] Real panel-copy transition structure

### Stage 3: bridge / stitching
- [x] Same / flip / abstain scaffold
- [x] Stitched block path
- [x] Final stitched site-level output
- [ ] Real bridge evidence from mature local blocks
- [ ] Better block construction than naive sign-run blockify

### Stage 4: evaluation
- [x] Toy scenarios
- [x] Run reports / TSV outputs
- [ ] HG002 chr20 truth benchmarking
- [ ] Calibration against truth
- [ ] Baseline comparison

## Current focus
Move from window-state scoring surrogate toward a true reduced-state local structured decoder.
