#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-readobs-json", required=True)
    ap.add_argument("--selected-reads-tsv", required=True)
    ap.add_argument("--haplotag-tsv", required=True)
    ap.add_argument("--superreads-tsv", required=True)
    ap.add_argument("--out-hybrid-json", required=True)
    ap.add_argument("--out-summary-json", required=True)
    ap.add_argument("--tagged-multiplier", type=float, default=1.20)
    ap.add_argument("--superread-base-weight", type=float, default=0.90)
    ap.add_argument("--superread-weight-step", type=float, default=0.08)
    ap.add_argument("--superread-weight-cap", type=float, default=1.35)
    args = ap.parse_args()

    py = sys.executable
    run([
        py, "python/awphase_py/read_frontend/build_hybrid_readset.py",
        "--raw-readobs-json", args.raw_readobs_json,
        "--selected-reads-tsv", args.selected_reads_tsv,
        "--haplotag-tsv", args.haplotag_tsv,
        "--superreads-tsv", args.superreads_tsv,
        "--out-hybrid-json", args.out_hybrid_json,
        "--out-summary-json", args.out_summary_json,
        "--tagged-multiplier", str(args.tagged_multiplier),
        "--superread-base-weight", str(args.superread_base_weight),
        "--superread-weight-step", str(args.superread_weight_step),
        "--superread-weight-cap", str(args.superread_weight_cap),
    ])

if __name__ == "__main__":
    main()
