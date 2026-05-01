#!/usr/bin/env python3
import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

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
        raise SystemExit(f"Could not auto-detect column. Candidates={candidates}; Columns={fieldnames}")
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

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def truthy(x):
    if x is None:
        return False
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}

def normalize_read_id(rid):
    s = str(rid).strip()
    # Remove common mate suffixes and whitespace tokens while keeping Illumina QNAME identity.
    s = s.split()[0]
    s = re.sub(r"[/._-]([12])$", "", s)
    return s

def parse_allele(raw):
    if raw is None:
        return 0
    s = str(raw).strip().lower()
    if s in {"1", "+1", "1.0", "alt", "alternate", "a1", "allele1", "hap1", "h1", "hp1"}:
        return 1
    if s in {"0", "-1", "-1.0", "ref", "reference", "a0", "allele0", "hap2", "h2", "hp2"}:
        return -1
    return 0

def load_phase6_calls(path):
    pos_to_comp = {}
    comp_to_positions = defaultdict(list)

    for r in read_tsv(path):
        pos = ival(r.get("pos"))
        st = ival(r.get("local_phase_state", r.get("phase_state", 0)))
        comp = str(r.get("phase6_component_id", r.get("block_id", ""))).strip()
        block = str(r.get("block_id", "")).strip()

        if comp == "":
            comp = block

        if pos <= 0 or st == 0:
            continue
        if comp == "" or comp.lower() in {"unassigned", "none", "na", "."}:
            continue

        pos_to_comp[pos] = {
            "component_id": comp,
            "state": st,
            "block_id": block,
        }
        comp_to_positions[comp].append(pos)

    for comp in comp_to_positions:
        comp_to_positions[comp].sort()

    return pos_to_comp, comp_to_positions

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--phase6-local-calls-tsv", required=True)
    ap.add_argument("--out-edges-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--chrom-start", type=int, default=None)
    ap.add_argument("--chrom-end", type=int, default=None)

    ap.add_argument("--read-id-col")
    ap.add_argument("--pos-col")
    ap.add_argument("--allele-col")
    ap.add_argument("--weight-col")

    ap.add_argument("--min-obs-weight", type=float, default=0.005)
    ap.add_argument("--max-obs-weight", type=float, default=10.0)
    ap.add_argument("--min-component-vote", type=float, default=0.02)
    ap.add_argument("--max-components-per-read", type=int, default=8)
    ap.add_argument("--include-ambiguous", action="store_true")
    args = ap.parse_args()

    pos_to_comp, comp_to_positions = load_phase6_calls(args.phase6_local_calls_tsv)

    if not pos_to_comp:
        raise SystemExit("No solved Phase 6 positions found.")

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

        # read_id -> component_id -> list of weighted orientation votes.
        # vote = phase6_state_at_site * observed_allele
        # This is the read's implied orientation relative to that component.
        per_read_comp_votes = defaultdict(lambda: defaultdict(list))

        rows_in = 0
        rows_used = 0
        skipped_ambiguous = 0
        skipped_region = 0
        skipped_unsolved = 0
        skipped_allele = 0
        skipped_weight = 0

        for r in reader:
            rows_in += 1

            if not args.include_ambiguous and any(truthy(r.get(c)) for c in ambiguous_cols):
                skipped_ambiguous += 1
                continue

            rid = normalize_read_id(r.get(read_col, ""))
            if not rid:
                continue

            pos = ival(r.get(pos_col), 0)
            if pos <= 0:
                continue

            if args.chrom_start is not None and pos < args.chrom_start:
                skipped_region += 1
                continue
            if args.chrom_end is not None and pos > args.chrom_end:
                skipped_region += 1
                continue

            pc = pos_to_comp.get(pos)
            if pc is None:
                skipped_unsolved += 1
                continue

            allele = parse_allele(r.get(allele_col))
            if allele == 0:
                skipped_allele += 1
                continue

            weight = 1.0
            if weight_col:
                weight = fval(r.get(weight_col), 1.0)
            if weight is None:
                weight = 1.0

            weight = abs(float(weight))
            if weight < args.min_obs_weight:
                skipped_weight += 1
                continue
            weight = min(weight, args.max_obs_weight)

            comp = pc["component_id"]
            phase_state = pc["state"]
            vote = phase_state * allele

            per_read_comp_votes[rid][comp].append((vote, weight, pos))
            rows_used += 1

    edge = defaultdict(lambda: {
        "plus_weight": 0.0,
        "minus_weight": 0.0,
        "plus_reads": 0,
        "minus_reads": 0,
        "read_examples": [],
    })

    reads_touching_components = 0
    reads_with_bridge = 0
    reads_skipped_too_many_components = 0

    for rid, comp_votes in per_read_comp_votes.items():
        comp_orient = {}

        for comp, votes in comp_votes.items():
            plus = sum(w for v, w, _ in votes if v == 1)
            minus = sum(w for v, w, _ in votes if v == -1)
            score = plus - minus

            if abs(score) < args.min_component_vote:
                continue

            sign = 1 if score > 0 else -1
            strength = abs(score)
            n_sites = len({p for _, _, p in votes})

            comp_orient[comp] = {
                "sign": sign,
                "strength": strength,
                "n_sites": n_sites,
            }

        if len(comp_orient) >= 1:
            reads_touching_components += 1

        comps = sorted(comp_orient)
        if len(comps) < 2:
            continue

        if len(comps) > args.max_components_per_read:
            reads_skipped_too_many_components += 1
            continue

        reads_with_bridge += 1

        for i in range(len(comps)):
            for j in range(i + 1, len(comps)):
                a = comps[i]
                b = comps[j]
                va = comp_orient[a]
                vb = comp_orient[b]

                # Required component relation:
                # orientation_b = orientation_a * edge_sign
                edge_sign = va["sign"] * vb["sign"]
                w = min(va["strength"], vb["strength"])

                key = (a, b)
                if edge_sign == 1:
                    edge[key]["plus_weight"] += w
                    edge[key]["plus_reads"] += 1
                else:
                    edge[key]["minus_weight"] += w
                    edge[key]["minus_reads"] += 1

                if len(edge[key]["read_examples"]) < 5:
                    edge[key]["read_examples"].append(rid)

    rows = []
    for (a, b), e in edge.items():
        plus = e["plus_weight"]
        minus = e["minus_weight"]
        total = plus + minus
        if total <= 0:
            continue

        sign = 1 if plus >= minus else -1
        support = max(plus, minus)
        conflict = min(plus, minus)
        margin = support - conflict
        reads_support = e["plus_reads"] if sign == 1 else e["minus_reads"]
        reads_conflict = e["minus_reads"] if sign == 1 else e["plus_reads"]

        rows.append({
            "component_a": a,
            "component_b": b,
            "edge_sign": sign,
            "support_weight": f"{support:.6f}",
            "conflict_weight": f"{conflict:.6f}",
            "margin_weight": f"{margin:.6f}",
            "total_weight": f"{total:.6f}",
            "support_reads": reads_support,
            "conflict_reads": reads_conflict,
            "plus_weight": f"{plus:.6f}",
            "minus_weight": f"{minus:.6f}",
            "plus_reads": e["plus_reads"],
            "minus_reads": e["minus_reads"],
            "read_examples": ",".join(e["read_examples"]),
        })

    rows.sort(
        key=lambda r: (
            float(r["margin_weight"]),
            int(r["support_reads"]),
            float(r["support_weight"]),
        ),
        reverse=True,
    )

    out = Path(args.out_edges_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "component_a", "component_b", "edge_sign",
        "support_weight", "conflict_weight", "margin_weight", "total_weight",
        "support_reads", "conflict_reads",
        "plus_weight", "minus_weight", "plus_reads", "minus_reads",
        "read_examples",
    ]

    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "obs_tsv": args.obs_tsv,
        "phase6_local_calls_tsv": args.phase6_local_calls_tsv,
        "out_edges_tsv": args.out_edges_tsv,
        "detected_columns": {
            "read_id": read_col,
            "pos": pos_col,
            "allele": allele_col,
            "weight": weight_col,
            "ambiguous_cols": ambiguous_cols,
        },
        "phase6_components": len(comp_to_positions),
        "phase6_solved_positions": len(pos_to_comp),
        "rows_in": rows_in,
        "rows_used": rows_used,
        "skipped_ambiguous": skipped_ambiguous,
        "skipped_region": skipped_region,
        "skipped_unsolved": skipped_unsolved,
        "skipped_allele": skipped_allele,
        "skipped_weight": skipped_weight,
        "reads_touching_components": reads_touching_components,
        "reads_with_bridge": reads_with_bridge,
        "reads_skipped_too_many_components": reads_skipped_too_many_components,
        "raw_edges": len(edge),
        "edges_written": len(rows),
        "min_obs_weight": args.min_obs_weight,
        "min_component_vote": args.min_component_vote,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
