#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import pysam

def load_variants(path, chrom=None, start=None, end=None):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    variants = {}
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue

        c = str(r.get("chrom", r.get("contig", chrom or "")))
        if chrom and c not in {"", chrom}:
            continue
        if start is not None and pos < start:
            continue
        if end is not None and pos > end:
            continue

        ref = str(r.get("ref_allele", r.get("ref", ""))).upper()
        alt = str(r.get("alt_allele", r.get("alt", ""))).upper()

        # Phase 6C v1 is SNP-only. Indels need realignment-aware handling.
        if len(ref) != 1 or len(alt) != 1:
            continue

        variants[pos] = {
            "pos": pos,
            "ref": ref,
            "alt": alt,
        }

    return variants

def base_weight(baseq, mapq, min_baseq, max_weight):
    # Conservative but simple. BaseQ controls most of the confidence; MAPQ downweights weak placements.
    b = max(0, baseq - min_baseq + 1)
    q_weight = min(1.0, b / 30.0)
    m_weight = min(1.0, max(0, mapq) / 60.0)
    return max(0.0, min(max_weight, q_weight * m_weight))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=15)
    ap.add_argument("--min-sites-per-fragment", type=int, default=2)
    ap.add_argument("--max-weight", type=float, default=1.0)
    ap.add_argument("--include-duplicates", action="store_true")
    ap.add_argument("--include-secondary", action="store_true")
    args = ap.parse_args()

    variants = load_variants(args.variant_json, args.chrom, args.start, args.end)
    if not variants:
        raise SystemExit("No SNP variants loaded from variant JSON.")

    bam = pysam.AlignmentFile(args.bam, "rb")
    if not bam.has_index():
        raise SystemExit(f"BAM/CRAM is not indexed: {args.bam}")

    # template_id -> pos -> list[(allele, weight)]
    template_obs = defaultdict(lambda: defaultdict(list))

    reads_in = 0
    reads_used = 0
    obs_seen = 0
    obs_used = 0
    skipped_mapq = 0
    skipped_flags = 0
    skipped_baseq = 0
    skipped_nonallelic_base = 0

    # pysam fetch is 0-based half-open.
    for read in bam.fetch(args.chrom, args.start - 1, args.end):
        reads_in += 1

        if read.is_unmapped:
            skipped_flags += 1
            continue
        if read.is_duplicate and not args.include_duplicates:
            skipped_flags += 1
            continue
        if (read.is_secondary or read.is_supplementary) and not args.include_secondary:
            skipped_flags += 1
            continue
        if read.mapping_quality < args.min_mapq:
            skipped_mapq += 1
            continue
        if read.query_sequence is None:
            continue

        qname = read.query_name
        seq = read.query_sequence.upper()
        quals = read.query_qualities

        read_had_obs = False

        # matches_only=False lets us ignore insertions/deletions cleanly.
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if qpos is None or rpos is None:
                continue

            pos = rpos + 1
            v = variants.get(pos)
            if v is None:
                continue

            obs_seen += 1

            if quals is not None:
                bq = int(quals[qpos])
            else:
                bq = 30

            if bq < args.min_baseq:
                skipped_baseq += 1
                continue

            base = seq[qpos]
            if base == v["alt"]:
                allele = 1
            elif base == v["ref"]:
                allele = -1
            else:
                skipped_nonallelic_base += 1
                continue

            w = base_weight(bq, read.mapping_quality, args.min_baseq, args.max_weight)
            if w <= 0:
                continue

            template_obs[qname][pos].append((allele, w))
            obs_used += 1
            read_had_obs = True

        if read_had_obs:
            reads_used += 1

    bam.close()

    fragments = []
    duplicate_sites_collapsed = 0
    site_votes_used = 0

    for qname, by_pos in template_obs.items():
        positions = []
        alleles = []
        weights = []

        for pos, votes in sorted(by_pos.items()):
            if len(votes) > 1:
                duplicate_sites_collapsed += 1

            plus = sum(w for a, w in votes if a == 1)
            minus = sum(w for a, w in votes if a == -1)

            if plus == minus:
                continue
            if plus > minus:
                allele = 1
                weight = plus - minus
            else:
                allele = -1
                weight = minus - plus

            if weight <= 0:
                continue

            positions.append(pos)
            alleles.append(allele)
            weights.append(round(min(args.max_weight, weight), 6))
            site_votes_used += 1

        if len(positions) >= args.min_sites_per_fragment:
            fragments.append({
                "fragment_id": qname,
                "n_sites": len(positions),
                "positions_json": json.dumps(positions),
                "alleles_json": json.dumps(alleles),
                "weights_json": json.dumps(weights),
                "fragment_score": round(sum(weights), 6),
            })

    fragments.sort(key=lambda r: (-int(r["n_sites"]), -float(r["fragment_score"]), r["fragment_id"]))

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = ["fragment_id", "n_sites", "positions_json", "alleles_json", "weights_json", "fragment_score"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(fragments)

    summary = {
        "bam": args.bam,
        "variant_json": args.variant_json,
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "snp_variants_loaded": len(variants),
        "reads_in": reads_in,
        "reads_used": reads_used,
        "obs_seen": obs_seen,
        "obs_used": obs_used,
        "skipped_flags": skipped_flags,
        "skipped_mapq": skipped_mapq,
        "skipped_baseq": skipped_baseq,
        "skipped_nonallelic_base": skipped_nonallelic_base,
        "templates_with_obs": len(template_obs),
        "fragments": len(fragments),
        "duplicate_sites_collapsed": duplicate_sites_collapsed,
        "site_votes_used": site_votes_used,
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "min_sites_per_fragment": args.min_sites_per_fragment,
        "out_tsv": args.out_tsv,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
