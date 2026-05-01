#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def ival(x, default=0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return int(float(x))
    except Exception:
        return default

def fval(x, default=0.0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return float(x)
    except Exception:
        return default

class ParityDSU:
    def __init__(self):
        self.parent = {}
        self.parity = {}  # parity to parent; 0 same, 1 flipped
        self.size = {}
        self.start = {}
        self.end = {}

    def add(self, x, start=None, end=None):
        if x not in self.parent:
            self.parent[x] = x
            self.parity[x] = 0
            self.size[x] = 1
            self.start[x] = start if start is not None else 10**18
            self.end[x] = end if end is not None else -1

    def find(self, x):
        self.add(x)
        if self.parent[x] == x:
            return x, 0
        root, p = self.find(self.parent[x])
        self.parity[x] ^= p
        self.parent[x] = root
        return self.parent[x], self.parity[x]

    def can_union(self, a, b, rel_parity):
        # rel_parity: 0 means same orientation, 1 means flipped.
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            return (pa ^ pb) == rel_parity
        return True

    def merged_span(self, a, b):
        ra, _ = self.find(a)
        rb, _ = self.find(b)
        if ra == rb:
            return self.end[ra] - self.start[ra] + 1, self.size[ra]
        s = min(self.start[ra], self.start[rb])
        e = max(self.end[ra], self.end[rb])
        n = self.size[ra] + self.size[rb]
        return e - s + 1, n

    def union(self, a, b, rel_parity):
        ra, pa = self.find(a)
        rb, pb = self.find(b)

        if ra == rb:
            return (pa ^ pb) == rel_parity

        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
            a, b = b, a
            pa, pb = pb, pa

        # Need pa ^ parity[rb] ^ pb == rel_parity
        self.parent[rb] = ra
        self.parity[rb] = pa ^ pb ^ rel_parity
        self.size[ra] += self.size[rb]
        self.start[ra] = min(self.start[ra], self.start[rb])
        self.end[ra] = max(self.end[ra], self.end[rb])
        return True

def block_stats_from_calls(rows):
    stats = {}
    for r in rows:
        block = str(r.get("block_id", "unassigned"))
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        pos = ival(r.get("pos"), None)
        if pos is None or state == 0 or block in {"", "0", "unassigned"}:
            continue
        if block not in stats:
            stats[block] = {"start": pos, "end": pos, "n": 0}
        stats[block]["start"] = min(stats[block]["start"], pos)
        stats[block]["end"] = max(stats[block]["end"], pos)
        stats[block]["n"] += 1
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--bridge-predictions-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--out-tsv", required=True)

    ap.add_argument("--min-prob", type=float, default=0.995)
    ap.add_argument("--min-confidence", type=float, default=0.50)
    ap.add_argument("--min-margin", type=float, default=50.0)
    ap.add_argument("--min-total-support", type=float, default=0.0)
    ap.add_argument("--min-effective-donors", type=float, default=0.0)
    ap.add_argument("--min-usable-anchors", type=int, default=2)
    ap.add_argument("--max-gap-bp", type=int, default=75000)

    ap.add_argument("--max-degree", type=int, default=1)
    ap.add_argument("--max-component-blocks", type=int, default=3)
    ap.add_argument("--max-component-span-bp", type=int, default=150000)
    ap.add_argument("--max-edges-per-window", type=int, default=999999)

    args = ap.parse_args()

    rows = read_tsv(args.local_calls_tsv)
    preds = read_tsv(args.bridge_predictions_tsv)

    stats = block_stats_from_calls(rows)

    dsu = ParityDSU()
    for b, st in stats.items():
        dsu.add(b, st["start"], st["end"])

    candidates = []
    for r in preds:
        if str(r.get("source_window")) != args.source_window:
            continue

        p = fval(r.get("xgb_oof_probability"), -1)
        conf = fval(r.get("bridge_confidence"), 0)
        margin = fval(r.get("bridge_margin"), 0)
        total = fval(r.get("bridge_total_support"), 0)
        eff = fval(r.get("bridge_effective_donors"), 0)
        usable_a = fval(r.get("bridge_a_usable_anchors"), 0)
        usable_b = fval(r.get("bridge_b_usable_anchors"), 0)
        gap = fval(r.get("gap_bp"), 10**18)

        if p < args.min_prob:
            continue
        if conf < args.min_confidence:
            continue
        if margin < args.min_margin:
            continue
        if total < args.min_total_support:
            continue
        if eff < args.min_effective_donors:
            continue
        if usable_a < args.min_usable_anchors or usable_b < args.min_usable_anchors:
            continue
        if gap > args.max_gap_bp:
            continue

        a = str(r.get("block_a", ""))
        b = str(r.get("block_b", ""))
        rel = ival(r.get("pred_relation", r.get("bridge_pred_relation", 0)), 0)

        if not a or not b or a == b or rel not in {-1, 1}:
            continue
        if a not in stats or b not in stats:
            continue

        rel_parity = 0 if rel == 1 else 1

        candidates.append({
            "a": a,
            "b": b,
            "rel": rel,
            "rel_parity": rel_parity,
            "prob": p,
            "conf": conf,
            "margin": margin,
            "total": total,
            "eff": eff,
            "gap": gap,
        })

    candidates.sort(
        key=lambda x: (
            x["prob"],
            x["conf"],
            x["margin"],
            x["total"],
            -x["gap"],
        ),
        reverse=True,
    )

    degree = defaultdict(int)
    accepted = []
    rejected_conflict = 0
    rejected_degree = 0
    rejected_span = 0
    rejected_component_size = 0

    for e in candidates:
        if len(accepted) >= args.max_edges_per_window:
            break

        a, b = e["a"], e["b"]

        if degree[a] >= args.max_degree or degree[b] >= args.max_degree:
            rejected_degree += 1
            continue

        if not dsu.can_union(a, b, e["rel_parity"]):
            rejected_conflict += 1
            continue

        span, nblocks = dsu.merged_span(a, b)
        if nblocks > args.max_component_blocks:
            rejected_component_size += 1
            continue
        if span > args.max_component_span_bp:
            rejected_span += 1
            continue

        ok = dsu.union(a, b, e["rel_parity"])
        if not ok:
            rejected_conflict += 1
            continue

        degree[a] += 1
        degree[b] += 1
        accepted.append(e)

    out_rows = []
    merged_sites = 0

    for r in rows:
        rr = dict(r)
        block = str(rr.get("block_id", "unassigned"))
        state = ival(rr.get("local_phase_state", rr.get("phase_state", 0)), 0)

        rr["phase9c_bridge_merged"] = 0
        rr["phase9c_original_block_id"] = block

        if block in stats and state != 0:
            root, parity = dsu.find(block)
            if root != block or dsu.size.get(root, 1) > 1:
                mult = -1 if parity == 1 else 1
                new_state = state * mult
                if "local_phase_state" in rr:
                    rr["local_phase_state"] = new_state
                if "phase_state" in rr:
                    rr["phase_state"] = new_state
                rr["block_id"] = "phase9c_forest_" + root
                rr["phase9c_bridge_merged"] = 1
                merged_sites += 1

        out_rows.append(rr)

    fields = list(rows[0].keys()) if rows else []
    for f in ["phase9c_bridge_merged", "phase9c_original_block_id"]:
        if f not in fields:
            fields.append(f)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    print({
        "source_window": args.source_window,
        "candidate_edges_after_filters": len(candidates),
        "accepted_edges": len(accepted),
        "merged_sites": merged_sites,
        "rejected_conflict": rejected_conflict,
        "rejected_degree": rejected_degree,
        "rejected_span": rejected_span,
        "rejected_component_size": rejected_component_size,
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
