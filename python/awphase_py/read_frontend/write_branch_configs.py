#!/usr/bin/env python3
import argparse
from pathlib import Path

def write_cfg(src_path, out_path, target_reads):
    txt = Path(src_path).read_text()
    lines = txt.splitlines()
    new_lines = []
    replaced = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("target_reads:"):
            prefix = line[:len(line) - len(line.lstrip())]
            new_lines.append(f"{prefix}target_reads: {target_reads}")
            replaced = True
        else:
            new_lines.append(line)
    if not replaced:
        raise SystemExit(f"Could not find target_reads in {src_path}")
    Path(out_path).write_text("\n".join(new_lines) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-config", required=True)
    ap.add_argument("--derived-stem", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    write_cfg(args.source_config, f"{args.out_prefix}.selected.yaml",   f"{args.derived_stem}_selected.json")
    write_cfg(args.source_config, f"{args.out_prefix}.tagged.yaml",     f"{args.derived_stem}_tagged.json")
    write_cfg(args.source_config, f"{args.out_prefix}.superreads.yaml", f"{args.derived_stem}_superreads.json")

if __name__ == "__main__":
    main()
