#!/usr/bin/env python3
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

def read_variant_json(path):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    variants = {}
    ordered = []
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        ref = str(r.get("ref_allele", r.get("ref", "")))
        alt = str(r.get("alt_allele", r.get("alt", "")))
        variants[pos] = {"pos": pos, "ref": ref, "alt": alt}
        ordered.append(pos)

    return variants, sorted(set(ordered))

def pick_col(fieldnames, user_value, candidates, required=True):
    if user_value:
        if user_value not in fieldnames:
            raise SystemExit(f"Requested column {user_value!r} not found. Columns={fieldnames}")
        return user_value

    lower = {c.lower(): c for c in fieldnames}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]

    if required:
        raise SystemExit(
            "Could not auto-detect required column. "
            f"Candidates={candidates}; Columns={fieldnames}"
        )
    return None

def fval(x, default=None):
    try:
        if x is None or str(x).strip() == "":
            return default
        v = float(x)
        if math.isnan(v):
            return default
        return v
    except Exception:
        return default

def parse_allele(raw, ref="", alt=""):
    if raw is None:
        return 0
    s = str(raw).strip()
    if s == "":
        return 0

    sl = s.lower()

    # Numeric encodings: 0=REF, 1=ALT, -1=REF/opposite.
    if sl in {"1", "+1", "1.0", "alt", "alternate", "a1", "allele1", "hap1", "h1", "hp1"}:
        return 1
    if sl in {"0", "-1", "-1.0", "ref", "reference", "a0", "allele0", "hap2", "h2", "hp2"}:
        return -1

    # Base/string comparison.
    if alt and s.upper() == alt.upper():
        return 1
    if ref and s.upper() == ref.upper():
        return -1

    return 0

def truthy(x):
    if x is None:
        return False
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--chrom-start", type=int, default=None)
    ap.add_argument("--chrom-end", type=int, default=None)

    ap.add_argument("--read-id-col")
    ap.add_argument("--pos-col")
    ap.add_argument("--allele-col")
    ap.add_argument("--weight-col")

    ap.add_argument("--min-weight", type=float, default=0.05)
    ap.add_argument("--max-weight", type=float, default=10.0)
    ap.add_argument("--min-sites-per-fragment", type=int, default=2)
    ap.add_argument("--min-site-margin", type=float, default=0.0)

    args = ap.parse_args()

    variants, variant_positions = read_variant_json(args.variant_json)
    variant_pos_set = set(variant_positions)

    with open(args.obs_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        read_col = pick_col(
            fields,
            args.read_id_col,
            ["read_id", "read_name", "qname", "query_name", "fragment_id", "name", "read"],
            required=True,
        )
        pos_col = pick_col(
            fields,
            args.pos_col,
            ["site_pos", "pos", "position", "variant_pos"],
            required=True,
        )
        allele_col = pick_col(
            fields,
            args.allele_col,
            ["obs_allele", "allele", "read_allele", "variant_allele", "allele_state", "state", "base"],
            required=True,
        )
        weight_col = pick_col(
            fields,
            args.weight_col,
            ["weight", "obs_weight", "allele_weight", "allele_confidence", "confidence", "conf", "score", "abs_score_delta", "delta"],
            required=False,
        )

        ambiguous_cols = [
            c for c in fields
            if c.lower() in {"ambiguous", "is_ambiguous", "low_confidence", "is_low_confidence", "filtered"}
        ]

        raw = defaultdict(lambda: defaultdict(list))

        rows_in = 0
        rows_kept = 0
        skipped_region = 0
        skipped_nonvariant = 0
        skipped_allele = 0
        skipped_weight = 0
        skipped_ambiguous = 0

        for r in reader:
            rows_in += 1

            if any(truthy(r.get(c)) for c in ambiguous_cols):
                skipped_ambiguous += 1
                continue

            rid = str(r.get(read_col, "")).strip()
            if not rid:
                continue

            try:
                pos = int(float(r.get(pos_col, "")))
            except Exception:
                continue

            if args.chrom_start is not None and pos < args.chrom_start:
                skipped_region += 1
                continue
            if args.chrom_end is not None and pos > args.chrom_end:
                skipped_region += 1
                continue
            if pos not in variant_pos_set:
                skipped_nonvariant += 1
                continue

            v = variants[pos]
            allele = parse_allele(r.get(allele_col), v.get("ref", ""), v.get("alt", ""))
            if allele == 0:
                skipped_allele += 1
                continue

            weight = 1.0
            if weight_col:
                weight = fval(r.get(weight_col), 1.0)

            if weight is None:
                weight = 1.0

            weight = abs(float(weight))
            if weight < args.min_weight:
                skipped_weight += 1
                continue
            weight = min(weight, args.max_weight)

            raw[rid][pos].append((allele, weight))
            rows_kept += 1

    fragments = []

    site_obs = 0
    site_used = 0
    duplicate_sites_collapsed = 0

    for rid, by_pos in raw.items():
        positions = []
        alleles = []
        weights = []

        for pos, obs in sorted(by_pos.items()):
            site_obs += len(obs)
            if len(obs) > 1:
                duplicate_sites_collapsed += 1

            plus = sum(w for a, w in obs if a == 1)
            minus = sum(w for a, w in obs if a == -1)
            margin = abs(plus - minus)

            if margin < args.min_site_margin:
                continue

            if plus > minus:
                allele = 1
                weight = plus - minus if margin > 0 else plus
            elif minus > plus:
                allele = -1
                weight = minus - plus if margin > 0 else minus
            else:
                continue

            weight = max(args.min_weight, min(args.max_weight, weight))

            positions.append(pos)
            alleles.append(allele)
            weights.append(round(weight, 6))
            site_used += 1

        if len(positions) >= args.min_sites_per_fragment:
            fragments.append({
                "fragment_id": rid,
                "n_sites": len(positions),
                "positions_json": json.dumps(positions),
                "alleles_json": json.dumps(alleles),
                "weights_json": json.dumps(weights),
                "fragment_score": round(sum(weights), 6),
            })

    fragments.sort(key=lambda r: (-int(r["n_sites"]), -float(r["fragment_score"]), r["fragment_id"]))

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="") as fh:
        fields = ["fragment_id", "n_sites", "positions_json", "alleles_json", "weights_json", "fragment_score"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(fragments)

    summary = {
        "obs_tsv": args.obs_tsv,
        "variant_json": args.variant_json,
        "out_tsv": args.out_tsv,
        "detected_columns": {
            "read_id": read_col,
            "pos": pos_col,
            "allele": allele_col,
            "weight": weight_col,
            "ambiguous_cols": ambiguous_cols,
        },
        "rows_in": rows_in,
        "rows_kept": rows_kept,
        "skipped_region": skipped_region,
        "skipped_nonvariant": skipped_nonvariant,
        "skipped_allele": skipped_allele,
        "skipped_weight": skipped_weight,
        "skipped_ambiguous": skipped_ambiguous,
        "raw_reads": len(raw),
        "fragments": len(fragments),
        "site_observations_seen_after_grouping": site_obs,
        "site_observations_used_after_grouping": site_used,
        "duplicate_read_sites_collapsed": duplicate_sites_collapsed,
        "min_sites_per_fragment": args.min_sites_per_fragment,
        "min_weight": args.min_weight,
        "max_weight": args.max_weight,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
