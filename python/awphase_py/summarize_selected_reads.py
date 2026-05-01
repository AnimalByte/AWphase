#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def as_list(obj, preferred_keys):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        for k in preferred_keys:
            if k in obj and isinstance(obj[k], list):
                return obj[k]
        for v in obj.values():
            if isinstance(v, list):
                return v
    raise ValueError("Could not find a list payload in JSON")

def variant_positions(variants_obj):
    variants = as_list(variants_obj, ["variants", "sites", "items", "records"])
    out = set()
    for v in variants:
        if isinstance(v, dict):
            if "pos" in v:
                out.add(int(v["pos"]))
            elif "site_pos" in v:
                out.add(int(v["site_pos"]))
    return out

def read_records(reads_obj):
    return as_list(reads_obj, ["reads", "observations", "items", "records"])

def dedup_best_observation_per_site(obs):
    best = {}
    for o in obs:
        site = int(o["site_pos"])
        new_score = int(o.get("baseq", 0)) + int(o.get("mapq", 0))
        prev = best.get(site)
        if prev is None:
            best[site] = o
        else:
            prev_score = int(prev.get("baseq", 0)) + int(prev.get("mapq", 0))
            if new_score > prev_score:
                best[site] = o
    return [best[k] for k in sorted(best)]

def bundle_reads(variants_set, reads):
    grouped = defaultdict(list)
    covered_before = set()

    for r in reads:
        if not isinstance(r, dict):
            continue
        if "site_pos" not in r or "read_id" not in r:
            continue
        site = int(r["site_pos"])
        allele = int(r.get("allele", 0))
        if site not in variants_set:
            continue
        if allele == 0:
            continue
        grouped[str(r["read_id"])].append(r)
        covered_before.add(site)

    bundles = []
    for read_id, obs in grouped.items():
        dedup = dedup_best_observation_per_site(obs)
        if not dedup:
            continue

        unique_sites = [int(x["site_pos"]) for x in dedup]
        informative_sites = len(unique_sites)
        multi_site_bonus = max(0, informative_sites - 1)

        qual_sum = 0.0
        for x in dedup:
            bq = min(max(float(x.get("baseq", 0)) / 40.0, 0.0), 1.0)
            mq = min(max(float(x.get("mapq", 0)) / 60.0, 0.0), 1.0)
            qual_sum += 0.5 + 0.5 * bq * mq

        score = 3.0 * multi_site_bonus + 1.0 * informative_sites + 0.25 * qual_sum

        bundles.append({
            "read_id": read_id,
            "observations": dedup,
            "unique_sites": unique_sites,
            "score": score,
            "informative_sites": informative_sites,
            "multi_site_bonus": multi_site_bonus,
        })

    bundles.sort(
        key=lambda b: (
            -int(b["multi_site_bonus"]),
            -int(b["informative_sites"]),
            -float(b["score"]),
            str(b["read_id"]),
        )
    )
    return bundles, len(covered_before)

def reconstruct_selection(variants_json, reads_json, max_reads, max_per_site, min_sites_per_read):
    variants_obj = load_json(variants_json)
    reads_obj = load_json(reads_json)

    vset = variant_positions(variants_obj)
    reads = read_records(reads_obj)

    bundles, covered_before = bundle_reads(vset, reads)

    site_cov = defaultdict(int)
    kept = []
    annotated = []

    for bundle in bundles:
        reason = "kept"
        selected = False

        if len(kept) >= max_reads:
            reason = "max_reads_reached"
        elif int(bundle["informative_sites"]) < min_sites_per_read:
            reason = "below_min_sites_per_read"
        else:
            contributes_undercovered = any(site_cov[p] < max_per_site for p in bundle["unique_sites"])
            if not contributes_undercovered:
                reason = "no_undercovered_site"
            else:
                selected = True
                for p in bundle["unique_sites"]:
                    if site_cov[p] < max_per_site:
                        site_cov[p] += 1
                kept.append(bundle)

        annotated.append({
            "read_id": bundle["read_id"],
            "selected": selected,
            "decision": reason,
            "informative_sites": bundle["informative_sites"],
            "unique_sites_count": len(bundle["unique_sites"]),
            "observations_count": len(bundle["observations"]),
            "multi_site_bonus": bundle["multi_site_bonus"],
            "score": round(bundle["score"], 6),
            "unique_sites": bundle["unique_sites"],
        })

    covered_after = sum(1 for _, c in site_cov.items() if c > 0)

    summary = {
        "input_observations": len(reads),
        "output_observations": sum(len(b["observations"]) for b in kept),
        "input_reads": len({str(r["read_id"]) for r in reads if isinstance(r, dict) and "read_id" in r}),
        "output_reads": len(kept),
        "target_sites": len(vset),
        "covered_sites_before": covered_before,
        "covered_sites_after": covered_after,
        "max_reads": max_reads,
        "max_per_site": max_per_site,
        "min_sites_per_read": min_sites_per_read,
    }

    return annotated, summary

def write_tsv(rows, out_path):
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "read_id",
                "selected",
                "decision",
                "informative_sites",
                "unique_sites_count",
                "observations_count",
                "multi_site_bonus",
                "score",
                "unique_sites",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for r in rows:
            row = dict(r)
            row["unique_sites"] = ",".join(map(str, row["unique_sites"]))
            writer.writerow(row)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", help="Existing selected_reads_debug.json, if you already have one")
    ap.add_argument("--reads-json", help="Read observations JSON")
    ap.add_argument("--variants-json", help="Variants JSON")
    ap.add_argument("--max-reads", type=int, default=20000)
    ap.add_argument("--max-per-site", type=int, default=12)
    ap.add_argument("--min-sites-per-read", type=int, default=3)
    ap.add_argument("--emit-json", help="Optional path to write reconstructed selected_reads_debug.json")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.json and Path(args.json).exists():
        obj = load_json(args.json)
        rows = as_list(obj, ["reads", "selected_reads", "items", "records"])
        write_tsv(rows, args.out)
        print(f"Wrote {len(rows)} rows to {args.out}")
        return

    if not args.reads_json or not args.variants_json:
        raise SystemExit(
            "selected_reads_debug.json does not exist yet. "
            "Use --reads-json and --variants-json to reconstruct it."
        )

    rows, summary = reconstruct_selection(
        variants_json=args.variants_json,
        reads_json=args.reads_json,
        max_reads=args.max_reads,
        max_per_site=args.max_per_site,
        min_sites_per_read=args.min_sites_per_read,
    )

    write_tsv(rows, args.out)

    if args.emit_json:
        payload = {"summary": summary, "reads": rows}
        out_json = Path(args.emit_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"Wrote reconstructed debug JSON to {out_json}")

    print(f"Wrote {len(rows)} rows to {args.out}")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
