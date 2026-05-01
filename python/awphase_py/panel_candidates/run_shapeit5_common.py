#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="phase_common")
    ap.add_argument("--input-bcf", required=True)
    ap.add_argument("--reference-bcf")
    ap.add_argument("--region", required=True)
    ap.add_argument("--map")
    ap.add_argument("--output-bcf", required=True)
    ap.add_argument("--thread", type=int, default=8)
    ap.add_argument("--log", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cmd = [
        args.binary,
        "--input", args.input_bcf,
        "--region", args.region,
        "--output", args.output_bcf,
        "--thread", str(args.thread),
    ]
    if args.reference_bcf:
        cmd.extend(["--reference", args.reference_bcf])
    if args.map:
        cmd.extend(["--map", args.map])

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    meta = {"tool": "shapeit5_phase_common", "cmd": cmd, "dry_run": args.dry_run}
    Path(args.log).write_text(json.dumps(meta, indent=2) + "\n")

    print(" ".join(cmd))
    if not args.dry_run:
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
