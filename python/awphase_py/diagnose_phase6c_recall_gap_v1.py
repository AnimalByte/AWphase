#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
import pysam
from collections import defaultdict

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
    return 0

def read_calls(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                pos = int(float(r["pos"]))
            except Exception:
                continue
            out[pos] = r
    return out

def load_variants(path):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj
    out = {}
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        ref = str(r.get("ref_allele", r.get("ref", ""))).upper()
        alt = str(r.get("alt_allele", r.get("alt", ""))).upper()
        vtype = "snp" if len(ref) == 1 and len(alt) == 1 else "indel_or_complex"
        out[pos] = {"pos": pos, "ref": ref, "alt": alt, "type": vtype}
    return out

def count_bam_observations(bam_path, chrom, positions):
    if not positions:
        return {}

    pos_set = set(positions)
    lo = min(pos_set)
    hi = max(pos_set)

    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = defaultdict(lambda: {"reads_covering": 0, "templates_covering": set()})

    for read in bam.fetch(chrom, lo - 1, hi):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if rpos is None:
                continue
            pos = rpos + 1
            if pos in pos_set:
                counts[pos]["reads_covering"] += 1
                counts[pos]["templates_covering"].add(qname)

    bam.close()

    return {
        p: {
            "reads_covering": counts[p]["reads_covering"],
            "templates_covering": len(counts[p]["templates_covering"]),
        }
        for p in pos_set
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6c-tsv", required=True)
    ap.add_argument("--whatshap-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    p6 = read_calls(args.phase6c_tsv)
    wh = read_calls(args.whatshap_tsv)
    variants = load_variants(args.variant_json)

    p6_nonzero = {p for p, r in p6.items() if state(r) != 0}
    wh_nonzero = {p for p, r in wh.items() if state(r) != 0}

    wh_only = sorted(wh_nonzero - p6_nonzero)
    p6_only = sorted(p6_nonzero - wh_nonzero)
    both = sorted(wh_nonzero & p6_nonzero)

    coverage = count_bam_observations(args.bam, args.chrom, wh_only)

    rows = []
    type_counts = defaultdict(int)

    for p in wh_only:
        v = variants.get(p, {"ref": "", "alt": "", "type": "missing_variant_json"})
        cov = coverage.get(p, {"reads_covering": 0, "templates_covering": 0})
        type_counts[v["type"]] += 1
        rows.append({
            "pos": p,
            "ref": v["ref"],
            "alt": v["alt"],
            "type": v["type"],
            "whatshap_state": state(wh[p]),
            "phase6c_state": state(p6.get(p, {})),
            "reads_covering": cov["reads_covering"],
            "templates_covering": cov["templates_covering"],
        })

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "pos", "ref", "alt", "type",
        "whatshap_state", "phase6c_state",
        "reads_covering", "templates_covering",
    ]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "phase6c_nonzero": len(p6_nonzero),
        "whatshap_nonzero": len(wh_nonzero),
        "both_nonzero": len(both),
        "whatshap_only_nonzero": len(wh_only),
        "phase6c_only_nonzero": len(p6_only),
        "whatshap_only_type_counts": dict(type_counts),
        "whatshap_only_with_any_read_coverage": sum(1 for r in rows if int(r["reads_covering"]) > 0),
        "whatshap_only_with_two_template_coverage": sum(1 for r in rows if int(r["templates_covering"]) >= 2),
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
