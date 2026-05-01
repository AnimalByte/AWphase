#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description="Run AWPhase's WhatsHap-inspired non-CLI subset pipeline")
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--max-reads", type=int, default=5000)
    ap.add_argument("--min-sites-per-read", type=int, default=2)
    ap.add_argument("--haplotag-min-sites", type=int, default=1)
    ap.add_argument("--haplotag-min-margin", type=float, default=5.0)
    ap.add_argument("--superread-min-support", type=int, default=2)
    args = ap.parse_args()

    wd = Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)

    run([
        sys.executable, "python/awphase_py/select_informative_reads.py",
        "--obs-tsv", args.obs_tsv,
        "--out-tsv", str(wd / "selected_reads.tsv"),
        "--out-summary", str(wd / "selected_reads.summary.json"),
        "--max-reads", str(args.max_reads),
        "--min-sites-per-read", str(args.min_sites_per_read),
    ])

    run([
        sys.executable, "python/awphase_py/haplotag_lite.py",
        "--local-calls-tsv", args.local_calls_tsv,
        "--obs-tsv", args.obs_tsv,
        "--selected-reads-tsv", str(wd / "selected_reads.tsv"),
        "--out-tsv", str(wd / "haplotag_lite.tsv"),
        "--out-summary", str(wd / "haplotag_lite.summary.json"),
        "--min-sites", str(args.haplotag_min_sites),
        "--min-margin", str(args.haplotag_min_margin),
    ])

    run([
        sys.executable, "python/awphase_py/build_superreads_from_haplotag.py",
        "--haplotag-tsv", str(wd / "haplotag_lite.tsv"),
        "--obs-tsv", args.obs_tsv,
        "--out-tsv", str(wd / "superreads.tsv"),
        "--out-summary", str(wd / "superreads.summary.json"),
        "--min-support", str(args.superread_min_support),
    ])


if __name__ == "__main__":
    main()
