#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", default="java")
    ap.add_argument("--jar", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ref")
    ap.add_argument("--map")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--xmx", default="32g")
    ap.add_argument("--log", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cmd = [
        args.java,
        f"-Xmx{args.xmx}",
        "-jar", args.jar,
        f"gt={args.gt}",
        f"out={args.out_prefix}",
        f"nthreads={args.threads}",
    ]
    if args.ref:
        cmd.append(f"ref={args.ref}")
    if args.map:
        cmd.append(f"map={args.map}")

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "tool": "beagle_phase",
        "cmd": cmd,
        "dry_run": args.dry_run,
    }
    Path(args.log).write_text(json.dumps(meta, indent=2) + "\n")

    print(" ".join(cmd))
    if not args.dry_run:
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
