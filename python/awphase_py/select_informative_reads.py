#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from statistics import median
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


def normalize_read_id(read_id: str) -> str:
    # WhatsHap treats paired-end as a single read with a hole in the middle.
    # In many BAMs query_name is already shared across mates, so this is usually a no-op.
    return read_id


def load_obs(path):
    by_read = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rid = normalize_read_id(row["read_id"])
            row["site_pos"] = as_int(row.get("site_pos"))
            row["allele"] = as_int(row.get("allele"))
            row["effective_baseq"] = as_int(row.get("effective_baseq"), as_int(row.get("baseq"), 0))
            row["mapq"] = as_int(row.get("mapq"), 0)
            row["weight_v2"] = as_float(row.get("weight_v2"), 0.0)
            row["weight_v3"] = as_float(row.get("weight_v3"), row["weight_v2"])
            row["is_ambiguous"] = as_int(row.get("is_ambiguous"), 0)
            by_read[rid].append(row)
    return by_read


def greedy_select(reads_by_id, max_reads, min_sites_per_read, favor_bridging=True):
    stats = []
    site_cov = defaultdict(int)

    for rid, obs in reads_by_id.items():
        obs = sorted(obs, key=lambda r: r["site_pos"])
        uniq_sites = sorted({r["site_pos"] for r in obs if r["allele"] in (-1, 1)})
        if len(uniq_sites) < min_sites_per_read:
            continue
        span = uniq_sites[-1] - uniq_sites[0] if len(uniq_sites) >= 2 else 0
        bridge_pairs = max(0, len(uniq_sites) - 1)
        ambig = sum(r["is_ambiguous"] for r in obs)
        mean_w = sum(r["weight_v2"] for r in obs) / max(1, len(obs))
        # WhatsHap-like read selection: prefer informative multi-site reads and long bridges,
        # but penalize already overserved sites during greedy selection.
        base_score = (len(uniq_sites) ** 2) + 0.002 * span + (2.0 * bridge_pairs if favor_bridging else 0.0)
        base_score += 0.25 * mean_w * len(obs)
        base_score -= 0.5 * ambig
        stats.append({
            "read_id": rid,
            "n_obs": len(obs),
            "n_sites": len(uniq_sites),
            "span_bp": span,
            "bridge_pairs": bridge_pairs,
            "ambiguous_obs": ambig,
            "mean_weight_v2": round(mean_w, 6),
            "base_score": base_score,
            "sites": uniq_sites,
        })

    selected = []
    remaining = stats[:]
    while remaining and len(selected) < max_reads:
        best = None
        best_score = None
        best_idx = None
        for i, rec in enumerate(remaining):
            novelty_penalty = 0.0
            for pos in rec["sites"]:
                novelty_penalty += min(3.0, 0.4 * site_cov[pos])
            score = rec["base_score"] - novelty_penalty
            if best is None or score > best_score:
                best = rec
                best_score = score
                best_idx = i
        picked = remaining.pop(best_idx)
        picked["greedy_score"] = round(best_score, 6)
        selected.append(picked)
        for pos in picked["sites"]:
            site_cov[pos] += 1

    summary = {
        "reads_considered": len(stats),
        "reads_selected": len(selected),
        "median_selected_sites": median([r["n_sites"] for r in selected]) if selected else 0,
        "median_selected_span_bp": median([r["span_bp"] for r in selected]) if selected else 0,
    }
    return selected, summary


def main():
    ap = argparse.ArgumentParser(description="WhatsHap-style informative read selection for AWPhase")
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--max-reads", type=int, default=5000)
    ap.add_argument("--min-sites-per-read", type=int, default=2)
    ap.add_argument("--no-favor-bridging", action="store_true")
    args = ap.parse_args()

    reads = load_obs(args.obs_tsv)
    selected, summary = greedy_select(
        reads_by_id=reads,
        max_reads=args.max_reads,
        min_sites_per_read=args.min_sites_per_read,
        favor_bridging=not args.no_favor_bridging,
    )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        fieldnames = ["read_id", "n_obs", "n_sites", "span_bp", "bridge_pairs", "ambiguous_obs", "mean_weight_v2", "base_score", "greedy_score"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in selected:
            writer.writerow({k: r[k] for k in fieldnames})

    Path(args.out_summary).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
