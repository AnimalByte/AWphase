#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
import pysam

def load_variants(path):
    obj = json.load(open(path))
    rows = obj.get("variants", obj) if isinstance(obj, dict) else obj

    out = []
    for r in rows:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        out.append({
            "pos": pos,
            "ref": str(r.get("ref_allele", r.get("ref", ""))).upper(),
            "alt": str(r.get("alt_allele", r.get("alt", ""))).upper(),
        })
    return out

def gt_to_state(gt, phased):
    if gt is None or len(gt) != 2 or None in gt:
        return 0

    # Only use phased heterozygous biallelic sites.
    if not phased:
        return 0

    if tuple(gt) == (0, 1):
        return 1
    if tuple(gt) == (1, 0):
        return -1
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phased-vcf", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--sample", default="HG002")
    args = ap.parse_args()

    variants = load_variants(args.variant_json)

    vf = pysam.VariantFile(args.phased_vcf)
    if args.sample not in vf.header.samples:
        raise SystemExit(f"Sample {args.sample} not found in {args.phased_vcf}. Samples={list(vf.header.samples)}")

    call_by_pos = {}

    for rec in vf.fetch():
        if rec.alts is None or len(rec.alts) != 1:
            continue

        pos = int(rec.pos)
        ref = str(rec.ref).upper()
        alt = str(rec.alts[0]).upper()

        sample = rec.samples[args.sample]
        gt = sample.get("GT")
        phased = bool(sample.phased)
        state = gt_to_state(gt, phased)

        ps = sample.get("PS")
        if ps is None:
            ps = pos if state != 0 else "unassigned"

        block_id = f"whatshap_{ps}" if state != 0 else "unassigned"

        call_by_pos[(pos, ref, alt)] = {
            "pos": pos,
            "block_id": block_id,
            "phase_state": state,
            "local_phase_state": state,
            "confidence": "1.000000" if state != 0 else "0.000000",
            "whatshap_gt": "/".join("." if x is None else str(x) for x in gt) if gt else "",
            "whatshap_phased": int(phased),
            "whatshap_ps": ps,
        }

    vf.close()

    out_rows = []
    for v in variants:
        key = (v["pos"], v["ref"], v["alt"])
        r = call_by_pos.get(key)
        if r is None:
            r = {
                "pos": v["pos"],
                "block_id": "unassigned",
                "phase_state": 0,
                "local_phase_state": 0,
                "confidence": "0.000000",
                "whatshap_gt": "",
                "whatshap_phased": 0,
                "whatshap_ps": "",
            }
        out_rows.append(r)

    fields = [
        "pos",
        "block_id",
        "phase_state",
        "local_phase_state",
        "confidence",
        "whatshap_gt",
        "whatshap_phased",
        "whatshap_ps",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    phased = sum(1 for r in out_rows if int(r["local_phase_state"]) != 0)
    print({"rows": len(out_rows), "phased": phased, "out": args.out_tsv})

if __name__ == "__main__":
    main()
