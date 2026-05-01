#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path
import sys

def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--derived-stem", required=True)
    ap.add_argument("--max-reads", type=int, default=5000)
    ap.add_argument("--min-sites-per-read", type=int, default=2)
    ap.add_argument("--min-haplotag-sites", type=int, default=1)
    ap.add_argument("--min-haplotag-margin", type=float, default=5.0)
    ap.add_argument("--min-superread-support", type=int, default=2)
    ap.add_argument("--reuse-readobs", action="store_true")
    args = ap.parse_args()

    repo = Path.cwd()
    workdir = repo / args.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    derived_stem = repo / args.derived_stem
    derived_stem.parent.mkdir(parents=True, exist_ok=True)

    readobs_json = workdir / "readobs_refaware_v3.json"
    readobs_summary = workdir / "readobs_refaware_v3.summary.json"
    site_tsv = workdir / "readobs_refaware_v3.site_debug.tsv"
    read_tsv = workdir / "readobs_refaware_v3.read_debug.tsv"
    obs_tsv = workdir / "readobs_refaware_v3.obs_debug.tsv"

    selected_tsv = workdir / "selected_reads.tsv"
    selected_summary = workdir / "selected_reads.summary.json"
    haplotag_tsv = workdir / "haplotag_lite.tsv"
    haplotag_summary = workdir / "haplotag_lite.summary.json"
    superreads_tsv = workdir / "superreads.tsv"
    superreads_summary = workdir / "superreads.summary.json"
    readsets_summary = workdir / "readsets.summary.json"

    selected_json = Path(str(derived_stem) + "_selected.json")
    tagged_json = Path(str(derived_stem) + "_tagged.json")
    superreads_json = Path(str(derived_stem) + "_superreads.json")

    py = sys.executable

    if not args.reuse_readobs or not readobs_json.exists():
        run([
            py, "python/awphase_py/build_readobs_from_bam_vcf_v2.py",
            "--bam", args.bam,
            "--vcf", args.vcf,
            "--fasta", args.fasta,
            "--region", args.region,
            "--sample", args.sample,
            "--out-json", str(readobs_json),
            "--out-summary", str(readobs_summary),
            "--out-site-tsv", str(site_tsv),
            "--out-read-tsv", str(read_tsv),
            "--out-obs-tsv", str(obs_tsv),
        ])

    run([
        py, "python/awphase_py/whatshap_core_subset_pipeline.py",
        "--local-calls-tsv", args.local_calls_tsv,
        "--obs-tsv", str(obs_tsv),
        "--workdir", str(workdir),
    ])

    run([
        py, "python/awphase_py/build_readsets_from_whatshap_subset.py",
        "--raw-readobs-json", str(readobs_json),
        "--selected-reads-tsv", str(selected_tsv),
        "--haplotag-tsv", str(haplotag_tsv),
        "--superreads-tsv", str(superreads_tsv),
        "--out-selected-json", str(selected_json),
        "--out-tagged-json", str(tagged_json),
        "--out-superreads-json", str(superreads_json),
        "--out-summary-json", str(readsets_summary),
    ])

    summary = {
        "workdir": str(workdir),
        "derived_outputs": {
            "selected": str(selected_json),
            "tagged": str(tagged_json),
            "superreads": str(superreads_json),
        }
    }
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
