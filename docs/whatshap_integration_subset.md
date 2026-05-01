# WhatsHap-inspired non-CLI subset added to AWPhase

This is **not** a literal port of the full WhatsHap core solver. WhatsHap's published/project structure mixes Python, C++, and Cython and includes a full MEC-style phasing engine, CLI, pedigree/polyploid features, and more.

This AWPhase subset ports the parts that cleanly map onto the existing system without rewriting the whole decoder:

1. **Informative read selection**
   - `python/awphase_py/select_informative_reads.py`
   - prefers reads covering multiple heterozygous sites and longer bridges

2. **Haplotag-lite**
   - `python/awphase_py/haplotag_lite.py`
   - assigns reads to `HP=1/2/UNK` and `PS=block_id` from AWPhase local calls

3. **Superread construction**
   - `python/awphase_py/build_superreads_from_haplotag.py`
   - builds consensus site observations per `(phase_set, haplotype)`

4. **Pipeline driver**
   - `python/awphase_py/whatshap_core_subset_pipeline.py`

## Example

```bash
python python/awphase_py/whatshap_core_subset_pipeline.py \
  --local-calls-tsv results/experiments/chr1_refaware_v2_on/local_calls.tsv \
  --obs-tsv results/experiments/chr1_refaware_v3_on/readobs_refaware_v3.obs_debug.tsv \
  --workdir results/experiments/chr1_whatshap_subset
```

## Expected outputs

- `selected_reads.tsv`
- `haplotag_lite.tsv`
- `superreads.tsv`
- corresponding `*.summary.json`

These are intended to feed the next AWPhase step: either superread-aware local decoding or block-boundary / reranking features.
