#!/usr/bin/env python3
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def as_int(v, default=0):
    try:
        return int(v)
    except Exception:
        try:
            return int(float(v))
        except Exception:
            return default


def as_float(v, default=0.0):
    try:
        return float(v)
    except Exception:
        return default


def load_tags(path):
    tags = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tags[row["read_id"]] = row
    return tags


def load_obs(path):
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row["site_pos"] = as_int(row.get("site_pos"))
            row["allele"] = as_int(row.get("allele"))
            row["weight_v2"] = as_float(row.get("weight_v2"), 0.0)
            row["effective_baseq"] = as_int(row.get("effective_baseq"), 0)
            row["is_ambiguous"] = as_int(row.get("is_ambiguous"), 0)
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description="Build WhatsHap-like superreads from haplotagged AWPhase reads")
    ap.add_argument("--haplotag-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--min-support", type=int, default=2)
    args = ap.parse_args()

    tags = load_tags(args.haplotag_tsv)
    obs = load_obs(args.obs_tsv)

    grouped = defaultdict(list)
    for r in obs:
        tag = tags.get(r["read_id"])
        if not tag:
            continue
        if tag["phase_set"] == "unassigned" or tag["hp"] not in ("1", "2"):
            continue
        if r["allele"] not in (-1, 1) or r["is_ambiguous"]:
            continue
        grouped[(tag["phase_set"], tag["hp"], r["site_pos"])].append(r)

    out_rows = []
    support_counts = []
    for (phase_set, hp, pos), rows in grouped.items():
        wt = Counter()
        qual_sum = Counter()
        for r in rows:
            allele = r["allele"]
            w = max(0.1, r["weight_v2"]) * max(1.0, r["effective_baseq"]) / 40.0
            wt[allele] += w
            qual_sum[allele] += r["effective_baseq"]
        if not wt:
            continue
        support = len(rows)
        if support < args.min_support:
            continue
        allele = 1 if wt[1] >= wt[-1] else -1
        ref_w = wt[-1]
        alt_w = wt[1]
        margin = abs(alt_w - ref_w)
        out_rows.append({
            "superread_id": f"{phase_set}_HP{hp}",
            "phase_set": phase_set,
            "hp": hp,
            "site_pos": pos,
            "allele": allele,
            "support_reads": support,
            "ref_weight": round(ref_w, 6),
            "alt_weight": round(alt_w, 6),
            "margin": round(margin, 6),
            "mean_baseq": round((qual_sum[1] + qual_sum[-1]) / max(1, support), 3),
        })
        support_counts.append(support)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        fieldnames = ["superread_id", "phase_set", "hp", "site_pos", "allele", "support_reads", "ref_weight", "alt_weight", "margin", "mean_baseq"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(out_rows, key=lambda r: (r["phase_set"], r["hp"], r["site_pos"])))

    summary = {
        "superread_rows": len(out_rows),
        "distinct_superreads": len({r["superread_id"] for r in out_rows}),
        "median_support_reads": sorted(support_counts)[len(support_counts)//2] if support_counts else 0,
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
