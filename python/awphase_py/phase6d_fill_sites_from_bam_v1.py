#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import pysam

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def fval(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def state(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def set_state(row, s):
    if "local_phase_state" in row:
        row["local_phase_state"] = int(s)
    if "phase_state" in row:
        row["phase_state"] = int(s)

def load_variants(path, chrom=None, start=None, end=None):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    variants = {}
    order = []

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

        # Phase6D v1 fills SNPs only. Indels need a local realignment-aware fill pass.
        if len(ref) != 1 or len(alt) != 1:
            continue

        variants[pos] = {"pos": pos, "ref": ref, "alt": alt}
        order.append(pos)

    return variants, sorted(set(order))

def read_local_calls(path):
    rows = []
    by_pos = {}

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        for r in reader:
            try:
                pos = int(float(r["pos"]))
            except Exception:
                continue
            rows.append(r)
            by_pos[pos] = r

    return rows, by_pos, fields

def comp_id(row):
    c = str(row.get("phase6_component_id", "")).strip()
    if not c:
        c = str(row.get("block_id", "")).strip()
    if c == "" or c.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return c

def base_weight(baseq, mapq, min_baseq, max_weight):
    # Keep this simple and conservative.
    b = max(0, baseq - min_baseq + 1)
    q_weight = min(1.0, b / 30.0)
    m_weight = min(1.0, max(0, mapq) / 60.0)
    return max(0.0, min(max_weight, q_weight * m_weight))

def collapse_template_obs(votes):
    # votes: pos -> list[(allele, weight)]
    out = {}
    for pos, xs in votes.items():
        plus = sum(w for a, w in xs if a == 1)
        minus = sum(w for a, w in xs if a == -1)
        if plus == minus:
            continue
        if plus > minus:
            out[pos] = (1, plus - minus)
        else:
            out[pos] = (-1, minus - plus)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--backbone-local-calls-tsv", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-fill-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-baseq", type=int, default=10)
    ap.add_argument("--max-weight", type=float, default=1.0)

    ap.add_argument("--min-anchor-sites", type=int, default=1)
    ap.add_argument("--min-anchor-margin-weight", type=float, default=0.10)
    ap.add_argument("--min-support-templates", type=int, default=2)
    ap.add_argument("--min-support-weight", type=float, default=0.50)
    ap.add_argument("--min-margin-weight", type=float, default=0.50)
    ap.add_argument("--min-confidence", type=float, default=0.95)
    ap.add_argument("--min-best-vs-second-margin", type=float, default=0.25)
    ap.add_argument("--max-fill-dist-bp", type=int, default=0)
    ap.add_argument("--allow-overwrite", action="store_true")
    args = ap.parse_args()

    variants, variant_order = load_variants(args.variant_json, args.chrom, args.start, args.end)
    base_rows, base_by_pos, base_fields = read_local_calls(args.backbone_local_calls_tsv)

    if not variants:
        raise SystemExit("No SNP variants loaded.")
    if not base_rows:
        raise SystemExit("No backbone local calls loaded.")

    # Existing solved backbone sites.
    phased = {}
    comp_span = defaultdict(lambda: [10**18, -1])
    comp_positions = defaultdict(list)

    for pos, r in base_by_pos.items():
        s = state(r)
        c = comp_id(r)
        if s == 0 or not c:
            continue
        phased[pos] = {"state": s, "component_id": c}
        comp_span[c][0] = min(comp_span[c][0], pos)
        comp_span[c][1] = max(comp_span[c][1], pos)
        comp_positions[c].append(pos)

    candidate_positions = set(variants)
    if not args.allow_overwrite:
        candidate_positions = {p for p in candidate_positions if p not in phased}

    bam = pysam.AlignmentFile(args.bam, "rb")
    if not bam.has_index():
        raise SystemExit(f"BAM/CRAM is not indexed: {args.bam}")

    # qname -> pos -> list[(allele, weight)]
    template_obs = defaultdict(lambda: defaultdict(list))

    reads_in = 0
    reads_used = 0
    obs_seen = 0
    obs_used = 0
    skipped_flags = 0
    skipped_mapq = 0
    skipped_baseq = 0
    skipped_nonallelic = 0

    for read in bam.fetch(args.chrom, args.start - 1, args.end):
        reads_in += 1

        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
            skipped_flags += 1
            continue
        if read.mapping_quality < args.min_mapq:
            skipped_mapq += 1
            continue
        if read.query_sequence is None:
            continue

        seq = read.query_sequence.upper()
        quals = read.query_qualities
        qname = read.query_name

        had = False

        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if qpos is None or rpos is None:
                continue

            pos = rpos + 1
            v = variants.get(pos)
            if v is None:
                continue

            obs_seen += 1

            bq = int(quals[qpos]) if quals is not None else 30
            if bq < args.min_baseq:
                skipped_baseq += 1
                continue

            base = seq[qpos]
            if base == v["alt"]:
                allele = 1
            elif base == v["ref"]:
                allele = -1
            else:
                skipped_nonallelic += 1
                continue

            w = base_weight(bq, read.mapping_quality, args.min_baseq, args.max_weight)
            if w <= 0:
                continue

            template_obs[qname][pos].append((allele, w))
            obs_used += 1
            had = True

        if had:
            reads_used += 1

    bam.close()

    # candidate_pos -> component -> support/conflict
    accum = defaultdict(lambda: defaultdict(lambda: {
        "plus_weight": 0.0,
        "minus_weight": 0.0,
        "plus_templates": set(),
        "minus_templates": set(),
        "anchor_templates": set(),
    }))

    templates_with_anchor = 0
    templates_with_candidate = 0
    templates_contributing_fill = 0

    for qname, votes in template_obs.items():
        obs = collapse_template_obs(votes)
        if not obs:
            continue

        # Determine orientation of this template relative to each solved component.
        comp_votes = defaultdict(list)
        candidates_here = []

        for pos, (allele, weight) in obs.items():
            if pos in phased:
                c = phased[pos]["component_id"]
                s = phased[pos]["state"]
                vote = s * allele
                comp_votes[c].append((vote, weight, pos))
            elif pos in candidate_positions:
                candidates_here.append((pos, allele, weight))

        if comp_votes:
            templates_with_anchor += 1
        if candidates_here:
            templates_with_candidate += 1
        if not comp_votes or not candidates_here:
            continue

        comp_orient = {}

        for c, xs in comp_votes.items():
            plus = sum(w for v, w, _ in xs if v == 1)
            minus = sum(w for v, w, _ in xs if v == -1)
            n_anchor_sites = len({p for _, _, p in xs})
            margin = abs(plus - minus)

            if n_anchor_sites < args.min_anchor_sites:
                continue
            if margin < args.min_anchor_margin_weight:
                continue

            orient = 1 if plus >= minus else -1
            comp_orient[c] = {
                "orient": orient,
                "anchor_margin": margin,
                "anchor_sites": n_anchor_sites,
            }

        if not comp_orient:
            continue

        contributed = False

        for pos, allele, weight in candidates_here:
            for c, info in comp_orient.items():
                lo, hi = comp_span[c]
                if args.max_fill_dist_bp >= 0:
                    if pos < lo - args.max_fill_dist_bp or pos > hi + args.max_fill_dist_bp:
                        continue

                pred_state = info["orient"] * allele
                a = accum[pos][c]
                if pred_state == 1:
                    a["plus_weight"] += weight
                    a["plus_templates"].add(qname)
                else:
                    a["minus_weight"] += weight
                    a["minus_templates"].add(qname)
                a["anchor_templates"].add(qname)
                contributed = True

        if contributed:
            templates_contributing_fill += 1

    fill_rows = []
    chosen = {}

    for pos, by_comp in accum.items():
        comp_candidates = []

        for c, a in by_comp.items():
            plus = a["plus_weight"]
            minus = a["minus_weight"]
            support = max(plus, minus)
            conflict = min(plus, minus)
            margin = support - conflict
            total = support + conflict
            conf = support / total if total > 0 else 0.0
            pred = 1 if plus >= minus else -1
            support_templates = len(a["plus_templates"] if pred == 1 else a["minus_templates"])
            conflict_templates = len(a["minus_templates"] if pred == 1 else a["plus_templates"])

            comp_candidates.append({
                "pos": pos,
                "component_id": c,
                "pred_state": pred,
                "support_weight": support,
                "conflict_weight": conflict,
                "margin_weight": margin,
                "confidence": conf,
                "support_templates": support_templates,
                "conflict_templates": conflict_templates,
                "anchor_templates": len(a["anchor_templates"]),
            })

        if not comp_candidates:
            continue

        comp_candidates.sort(
            key=lambda x: (
                x["margin_weight"],
                x["support_templates"],
                x["support_weight"],
                -x["conflict_weight"],
            ),
            reverse=True,
        )

        best = comp_candidates[0]
        second_margin = comp_candidates[1]["margin_weight"] if len(comp_candidates) > 1 else 0.0
        best_vs_second = best["margin_weight"] - second_margin

        best["best_vs_second_margin"] = best_vs_second
        best["n_component_candidates"] = len(comp_candidates)

        accepted = (
            best["support_templates"] >= args.min_support_templates
            and best["support_weight"] >= args.min_support_weight
            and best["margin_weight"] >= args.min_margin_weight
            and best["confidence"] >= args.min_confidence
            and best_vs_second >= args.min_best_vs_second_margin
        )

        best["accepted"] = int(accepted)
        fill_rows.append(best)

        if accepted:
            chosen[pos] = best

    # Apply fills to output rows.
    out_rows = []
    filled = 0

    for r in base_rows:
        rr = dict(r)
        pos = ival(rr.get("pos"))

        if pos in chosen and (args.allow_overwrite or state(rr) == 0):
            f = chosen[pos]
            c = f["component_id"]

            set_state(rr, f["pred_state"])
            rr["block_id"] = c
            rr["phase6d_filled"] = 1
            rr["phase6d_fill_component_id"] = c
            rr["phase6d_fill_confidence"] = f"{f['confidence']:.6f}"
            rr["phase6d_fill_support_weight"] = f"{f['support_weight']:.6f}"
            rr["phase6d_fill_conflict_weight"] = f"{f['conflict_weight']:.6f}"
            rr["phase6d_fill_margin_weight"] = f"{f['margin_weight']:.6f}"
            rr["phase6d_fill_support_templates"] = f["support_templates"]
            filled += 1
        else:
            rr.setdefault("phase6d_filled", 0)
            rr.setdefault("phase6d_fill_component_id", "")
            rr.setdefault("phase6d_fill_confidence", "")
            rr.setdefault("phase6d_fill_support_weight", "")
            rr.setdefault("phase6d_fill_conflict_weight", "")
            rr.setdefault("phase6d_fill_margin_weight", "")
            rr.setdefault("phase6d_fill_support_templates", "")

        out_rows.append(rr)

    # Write local calls.
    fields = []
    seen = set()
    for r in out_rows:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # Write fill candidates.
    fill_fields = [
        "pos", "component_id", "pred_state", "accepted",
        "support_weight", "conflict_weight", "margin_weight", "confidence",
        "support_templates", "conflict_templates", "anchor_templates",
        "best_vs_second_margin", "n_component_candidates",
    ]

    with open(args.out_fill_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fill_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(fill_rows)

    summary = {
        "bam": args.bam,
        "backbone_local_calls_tsv": args.backbone_local_calls_tsv,
        "variant_json": args.variant_json,
        "chrom": args.chrom,
        "start": args.start,
        "end": args.end,
        "variants_loaded_snp": len(variants),
        "backbone_phased_sites": len(phased),
        "candidate_positions": len(candidate_positions),
        "components": len(comp_positions),
        "reads_in": reads_in,
        "reads_used": reads_used,
        "obs_seen": obs_seen,
        "obs_used": obs_used,
        "skipped_flags": skipped_flags,
        "skipped_mapq": skipped_mapq,
        "skipped_baseq": skipped_baseq,
        "skipped_nonallelic": skipped_nonallelic,
        "templates_with_obs": len(template_obs),
        "templates_with_anchor": templates_with_anchor,
        "templates_with_candidate": templates_with_candidate,
        "templates_contributing_fill": templates_contributing_fill,
        "fill_candidate_sites": len(fill_rows),
        "filled_sites": filled,
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "min_anchor_sites": args.min_anchor_sites,
        "min_anchor_margin_weight": args.min_anchor_margin_weight,
        "min_support_templates": args.min_support_templates,
        "min_support_weight": args.min_support_weight,
        "min_margin_weight": args.min_margin_weight,
        "min_confidence": args.min_confidence,
        "min_best_vs_second_margin": args.min_best_vs_second_margin,
        "max_fill_dist_bp": args.max_fill_dist_bp,
        "allow_overwrite": args.allow_overwrite,
        "out_tsv": args.out_tsv,
        "out_fill_tsv": args.out_fill_tsv,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
