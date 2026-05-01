#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    obj = json.load(open(args.in_json))
    rows = obj.get("variants", obj) if isinstance(obj, dict) else obj

    out = []
    for r in rows:
        try:
            pos = int(r["pos"])
        except Exception:
            continue

        c = str(r.get("chrom", r.get("contig", args.chrom)))
        if c and c != args.chrom:
            continue
        if args.start <= pos <= args.end:
            out.append(r)

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2) + "\n")

    print({
        "in_json": args.in_json,
        "out_json": args.out_json,
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "rows": len(out),
    })

if __name__ == "__main__":
    main()
