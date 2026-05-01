#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import defaultdict

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                fields.append(k)
                seen.add(k)

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def fval(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def state(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def set_state(row, s):
    if "local_phase_state" in row:
        row["local_phase_state"] = int(s)
    if "phase_state" in row:
        row["phase_state"] = int(s)

def comp_id(row):
    c = str(row.get("phase6_component_id", "")).strip()
    if c == "":
        c = str(row.get("block_id", "")).strip()
    if c == "" or c.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return c

class SignedDSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}
        self.sign_to_parent = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            self.sign_to_parent[x] = 1

    def find(self, x):
        self.add(x)
        if self.parent[x] == x:
            return x, 1

        root, s_parent = self.find(self.parent[x])
        s = self.sign_to_parent[x] * s_parent
        self.parent[x] = root
        self.sign_to_parent[x] = s
        return root, s

    def union(self, a, b, relation):
        # relation means orient_b = orient_a * relation
        self.add(a)
        self.add(b)

        ra, sa = self.find(a)  # orient_a = sa * orient_ra
        rb, sb = self.find(b)  # orient_b = sb * orient_rb

        if ra == rb:
            current = sa * sb
            return current == relation

        # Attach lower rank under higher rank.
        # Need sign_to_parent[rb] if parent[rb] = ra:
        # sb * orient_rb = sa * orient_ra * relation
        # orient_rb = sa * relation * sb * orient_ra
        if self.rank[ra] < self.rank[rb]:
            # Attach ra under rb instead.
            # orient_a = sa * orient_ra
            # orient_b = sb * orient_rb
            # orient_b = orient_a * relation
            # sb * orient_rb = sa * orient_ra * relation
            # orient_ra = sb * relation * sa * orient_rb
            self.parent[ra] = rb
            self.sign_to_parent[ra] = sb * relation * sa
        else:
            self.parent[rb] = ra
            self.sign_to_parent[rb] = sa * relation * sb
            if self.rank[ra] == self.rank[rb]:
                self.rank[ra] += 1

        return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6-local-calls-tsv", required=True)
    ap.add_argument("--bridge-edges-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-used-edges-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--min-support-reads", type=int, default=2)
    ap.add_argument("--min-support-weight", type=float, default=0.25)
    ap.add_argument("--min-margin-weight", type=float, default=0.10)
    ap.add_argument("--max-conflict-weight", type=float, default=999999.0)
    ap.add_argument("--max-conflict-read-frac", type=float, default=0.50)
    args = ap.parse_args()

    calls = read_tsv(args.phase6_local_calls_tsv)
    edges = read_tsv(args.bridge_edges_tsv)

    components = set()
    comp_positions = defaultdict(list)

    for r in calls:
        c = comp_id(r)
        s = state(r)
        if c and s != 0:
            components.add(c)
            comp_positions[c].append(ival(r.get("pos")))

    dsu = SignedDSU()
    for c in components:
        dsu.add(c)

    filtered_edges = []

    for e in edges:
        a = e["component_a"]
        b = e["component_b"]
        sign = ival(e.get("edge_sign"), 0)

        support_reads = ival(e.get("support_reads"), 0)
        conflict_reads = ival(e.get("conflict_reads"), 0)
        support_weight = fval(e.get("support_weight"), 0.0)
        conflict_weight = fval(e.get("conflict_weight"), 0.0)
        margin_weight = fval(e.get("margin_weight"), 0.0)

        if sign not in {-1, 1}:
            continue
        if a not in components or b not in components:
            continue
        if support_reads < args.min_support_reads:
            continue
        if support_weight < args.min_support_weight:
            continue
        if margin_weight < args.min_margin_weight:
            continue
        if conflict_weight > args.max_conflict_weight:
            continue

        denom = support_reads + conflict_reads
        conflict_frac = conflict_reads / denom if denom else 0.0
        if conflict_frac > args.max_conflict_read_frac:
            continue

        filtered_edges.append({
            "component_a": a,
            "component_b": b,
            "edge_sign": sign,
            "support_reads": support_reads,
            "conflict_reads": conflict_reads,
            "support_weight": support_weight,
            "conflict_weight": conflict_weight,
            "margin_weight": margin_weight,
            "conflict_read_frac": conflict_frac,
            "raw": e,
        })

    filtered_edges.sort(
        key=lambda e: (
            e["margin_weight"],
            e["support_reads"],
            e["support_weight"],
            -e["conflict_weight"],
        ),
        reverse=True,
    )

    used_rows = []
    used_edges = 0
    conflict_edges = 0

    for e in filtered_edges:
        ok = dsu.union(e["component_a"], e["component_b"], e["edge_sign"])
        status = "used" if ok else "conflict_rejected"
        if ok:
            used_edges += 1
        else:
            conflict_edges += 1

        row = dict(e["raw"])
        row["phase6b_status"] = status
        row["phase6b_conflict_read_frac"] = f"{e['conflict_read_frac']:.6f}"
        used_rows.append(row)

    # Root/orientation assignment.
    root_to_index = {}
    comp_to_root = {}
    comp_to_orient = {}

    for c in sorted(components, key=lambda x: min(comp_positions[x]) if comp_positions[x] else 10**18):
        root, orient = dsu.find(c)
        comp_to_root[c] = root
        comp_to_orient[c] = orient
        if root not in root_to_index:
            root_to_index[root] = len(root_to_index) + 1

    out = []

    for r in calls:
        rr = dict(r)
        c = comp_id(rr)
        s = state(rr)

        if c and s != 0:
            root = comp_to_root.get(c, c)
            orient = comp_to_orient.get(c, 1)
            new_state = s * orient
            block = f"phase6b_block_{root_to_index[root]}"

            set_state(rr, new_state)
            rr["block_id"] = block
            rr["phase6b_original_component_id"] = c
            rr["phase6b_root_component_id"] = root
            rr["phase6b_component_orientation"] = orient
        else:
            rr["phase6b_original_component_id"] = c
            rr["phase6b_root_component_id"] = ""
            rr["phase6b_component_orientation"] = ""

        out.append(rr)

    write_tsv(args.out_tsv, out)
    write_tsv(args.out_used_edges_tsv, used_rows)

    root_members = defaultdict(list)
    for c, root in comp_to_root.items():
        root_members[root].append(c)

    stitched_blocks = len(root_members)
    multi_component_blocks = sum(1 for comps in root_members.values() if len(comps) > 1)
    largest_component_merge = max((len(comps) for comps in root_members.values()), default=0)

    summary = {
        "phase6_local_calls_tsv": args.phase6_local_calls_tsv,
        "bridge_edges_tsv": args.bridge_edges_tsv,
        "components": len(components),
        "raw_edges": len(edges),
        "filtered_edges": len(filtered_edges),
        "used_edges": used_edges,
        "conflict_edges": conflict_edges,
        "stitched_blocks": stitched_blocks,
        "multi_component_blocks": multi_component_blocks,
        "largest_component_merge": largest_component_merge,
        "min_support_reads": args.min_support_reads,
        "min_support_weight": args.min_support_weight,
        "min_margin_weight": args.min_margin_weight,
        "max_conflict_weight": args.max_conflict_weight,
        "max_conflict_read_frac": args.max_conflict_read_frac,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
