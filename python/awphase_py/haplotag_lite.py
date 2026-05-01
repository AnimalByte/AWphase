#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
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


def load_local_calls(path):
    blocks = defaultdict(dict)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = as_int(row.get("pos"))
            block_id = row.get("block_id", "unassigned")
            state = as_int(row.get("local_phase_state", row.get("phase_state", 0)))
            conf = as_float(row.get("confidence", 1.0), 1.0)
            if pos <= 0 or state == 0 or block_id == "unassigned":
                continue
            blocks[block_id][pos] = {"state": state, "confidence": conf}
    return blocks


def load_obs(path, keep_selected=None):
    by_read = defaultdict(list)
    keep = None
    if keep_selected:
        keep = set()
        with open(keep_selected) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                keep.add(row["read_id"])
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rid = row["read_id"]
            if keep is not None and rid not in keep:
                continue
            row["site_pos"] = as_int(row.get("site_pos"))
            row["allele"] = as_int(row.get("allele"))
            row["weight_v2"] = as_float(row.get("weight_v2"), 0.0)
            row["is_ambiguous"] = as_int(row.get("is_ambiguous"), 0)
            row["effective_baseq"] = as_int(row.get("effective_baseq"), 0)
            by_read[rid].append(row)
    return by_read


def haplotag_read(obs, blocks, cost_scale="baseq", min_sites=1, min_margin=1.0):
    best = None
    for block_id, phase_map in blocks.items():
        cost_hp1 = 0.0
        cost_hp2 = 0.0
        used = 0
        used_pos = []
        for r in obs:
            pos = r["site_pos"]
            if pos not in phase_map:
                continue
            if r["allele"] not in (-1, 1):
                continue
            if r["is_ambiguous"]:
                continue
            state = phase_map[pos]["state"]
            conf = phase_map[pos]["confidence"]
            # AWPhase state convention: treat state=+1 as ALT on HP1 and REF on HP2.
            hp1_expected = state
            hp2_expected = -state
            if cost_scale == "baseq":
                w = max(1.0, float(r["effective_baseq"]))
            else:
                w = max(0.25, float(r["weight_v2"])) * 40.0
            w *= max(0.25, conf)
            if r["allele"] != hp1_expected:
                cost_hp1 += w
            if r["allele"] != hp2_expected:
                cost_hp2 += w
            used += 1
            used_pos.append(pos)
        if used < min_sites:
            continue
        margin = abs(cost_hp1 - cost_hp2)
        hp = "UNK"
        if margin >= min_margin:
            hp = "1" if cost_hp1 < cost_hp2 else "2"
        rec = {
            "phase_set": block_id,
            "hp": hp,
            "sites_used": used,
            "cost_hp1": round(cost_hp1, 6),
            "cost_hp2": round(cost_hp2, 6),
            "margin": round(margin, 6),
            "used_positions": used_pos,
        }
        if best is None or rec["margin"] > best["margin"]:
            best = rec
    if best is None:
        return {
            "phase_set": "unassigned",
            "hp": "UNK",
            "sites_used": 0,
            "cost_hp1": 0.0,
            "cost_hp2": 0.0,
            "margin": 0.0,
            "used_positions": [],
        }
    return best


def main():
    ap = argparse.ArgumentParser(description="WhatsHap-inspired haplotagging for AWPhase read observations")
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--selected-reads-tsv")
    ap.add_argument("--min-sites", type=int, default=1)
    ap.add_argument("--min-margin", type=float, default=5.0)
    ap.add_argument("--cost-scale", choices=["baseq", "weight_v2"], default="baseq")
    args = ap.parse_args()

    blocks = load_local_calls(args.local_calls_tsv)
    by_read = load_obs(args.obs_tsv, keep_selected=args.selected_reads_tsv)

    rows = []
    summary = defaultdict(int)
    for rid, obs in by_read.items():
        rec = haplotag_read(obs, blocks, cost_scale=args.cost_scale, min_sites=args.min_sites, min_margin=args.min_margin)
        rows.append({
            "read_id": rid,
            "phase_set": rec["phase_set"],
            "hp": rec["hp"],
            "sites_used": rec["sites_used"],
            "cost_hp1": rec["cost_hp1"],
            "cost_hp2": rec["cost_hp2"],
            "margin": rec["margin"],
        })
        summary[f"hp_{rec['hp']}"] += 1
        if rec["phase_set"] != "unassigned":
            summary["tagged_reads"] += 1
    summary["total_reads"] = len(rows)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        fieldnames = ["read_id", "phase_set", "hp", "sites_used", "cost_hp1", "cost_hp2", "margin"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    Path(args.out_summary).write_text(json.dumps(dict(summary), indent=2))
    print(json.dumps(dict(summary), indent=2))


if __name__ == "__main__":
    main()
