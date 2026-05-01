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
        self.parity = {}
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
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            return (pa ^ pb) == rel_parity
        return True

    def merged_stats(self, a, b):
        ra, _ = self.find(a)
        rb, _ = self.find(b)

        if ra == rb:
            span = self.end[ra] - self.start[ra] + 1
            return span, self.size[ra], 0

        old_max = max(
            self.end[ra] - self.start[ra] + 1,
            self.end[rb] - self.start[rb] + 1,
        )
        new_span = max(self.end[ra], self.end[rb]) - min(self.start[ra], self.start[rb]) + 1
        gain = new_span - old_max
        nblocks = self.size[ra] + self.size[rb]
        return new_span, nblocks, gain

    def union(self, a, b, rel_parity):
        ra, pa = self.find(a)
        rb, pb = self.find(b)

        if ra == rb:
            return (pa ^ pb) == rel_parity

        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
            a, b = b, a
            pa, pb = pb, pa

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


def relation_from_row(r):
    rel = ival(r.get("pred_relation", r.get("bridge_pred_relation", 0)), 0)
    if rel not in {-1, 1}:
        rel = ival(r.get("bridge_pred_relation", 0), 0)
    return rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--bridge-predictions-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--out-tsv", required=True)

    ap.add_argument("--min-prob", type=float, default=0.995)
    ap.add_argument("--min-confidence", type=float, default=0.50)
    ap.add_argument("--min-margin", type=float, default=50.0)
    ap.add_argument("--min-span-gain-bp", type=int, default=10000)
    ap.add_argument("--min-merged-span-bp", type=int, default=25000)
    ap.add_argument("--min-total-support", type=float, default=0.0)
    ap.add_argument("--min-effective-donors", type=float, default=0.0)
    ap.add_argument("--min-usable-anchors", type=int, default=2)
    ap.add_argument("--max-gap-bp", type=int, default=500000)
    ap.add_argument("--max-degree", type=int, default=1)
    ap.add_argument("--max-component-blocks", type=int, default=6)
    ap.add_argument("--max-component-span-bp", type=int, default=1000000)
    ap.add_argument("--max-edges-per-window", type=int, default=999999)
    ap.add_argument("--sort-mode", choices=["prob_first", "utility_first", "hybrid"], default="hybrid")
    ap.add_argument("--utility-weight", type=float, default=0.00001)

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

        p = fval(r.get("xgb_oof_probability"), -1.0)
        conf = fval(r.get("bridge_confidence"), 0.0)
        margin = fval(r.get("bridge_margin"), 0.0)
        total = fval(r.get("bridge_total_support"), 0.0)
        eff = fval(r.get("bridge_effective_donors"), 0.0)
        usable_a = fval(r.get("bridge_a_usable_anchors"), 0.0)
        usable_b = fval(r.get("bridge_b_usable_anchors"), 0.0)
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
        rel = relation_from_row(r)

        if not a or not b or a == b or rel not in {-1, 1}:
            continue
        if a not in stats or b not in stats:
            continue

        span_a = stats[a]["end"] - stats[a]["start"] + 1
        span_b = stats[b]["end"] - stats[b]["start"] + 1
        merged_span = max(stats[a]["end"], stats[b]["end"]) - min(stats[a]["start"], stats[b]["start"]) + 1
        span_gain = merged_span - max(span_a, span_b)

        if span_gain < args.min_span_gain_bp:
            continue
        if merged_span < args.min_merged_span_bp:
            continue

        rel_parity = 0 if rel == 1 else 1
        utility_score = span_gain * args.utility_weight

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
            "span_gain": span_gain,
            "merged_span": merged_span,
            "utility_score": utility_score,
        })

    if args.sort_mode == "prob_first":
        candidates.sort(key=lambda x: (x["prob"], x["conf"], x["margin"], x["span_gain"], -x["gap"]), reverse=True)
    elif args.sort_mode == "utility_first":
        candidates.sort(key=lambda x: (x["span_gain"], x["prob"], x["conf"], x["margin"], -x["gap"]), reverse=True)
    else:
        candidates.sort(
            key=lambda x: (
                x["prob"] + x["utility_score"],
                x["conf"],
                x["margin"],
                x["span_gain"],
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
    rejected_dynamic_gain = 0

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

        dyn_span, dyn_blocks, dyn_gain = dsu.merged_stats(a, b)

        if dyn_blocks > args.max_component_blocks:
            rejected_component_size += 1
            continue

        if dyn_span > args.max_component_span_bp:
            rejected_span += 1
            continue

        if dyn_gain < args.min_span_gain_bp:
            rejected_dynamic_gain += 1
            continue

        ok = dsu.union(a, b, e["rel_parity"])
        if not ok:
            rejected_conflict += 1
            continue

        degree[a] += 1
        degree[b] += 1
        e["dynamic_span"] = dyn_span
        e["dynamic_gain"] = dyn_gain
        accepted.append(e)

    out_rows = []
    merged_sites = 0

    for r in rows:
        rr = dict(r)
        block = str(rr.get("block_id", "unassigned"))
        state = ival(rr.get("local_phase_state", rr.get("phase_state", 0)), 0)

        rr["phase9d_bridge_merged"] = 0
        rr["phase9d_original_block_id"] = block

        if block in stats and state != 0:
            root, parity = dsu.find(block)

            if root != block or dsu.size.get(root, 1) > 1:
                mult = -1 if parity == 1 else 1
                new_state = state * mult

                if "local_phase_state" in rr:
                    rr["local_phase_state"] = new_state
                if "phase_state" in rr:
                    rr["phase_state"] = new_state

                rr["block_id"] = "phase9d_forest_" + root
                rr["phase9d_bridge_merged"] = 1
                merged_sites += 1

        out_rows.append(rr)

    fields = list(rows[0].keys()) if rows else []
    for f in ["phase9d_bridge_merged", "phase9d_original_block_id"]:
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
        "accepted_span_gain_sum": int(sum(e.get("dynamic_gain", 0) for e in accepted)),
        "accepted_mean_span_gain": (
            sum(e.get("dynamic_gain", 0) for e in accepted) / len(accepted)
        ) if accepted else 0,
        "rejected_conflict": rejected_conflict,
        "rejected_degree": rejected_degree,
        "rejected_span": rejected_span,
        "rejected_component_size": rejected_component_size,
        "rejected_dynamic_gain": rejected_dynamic_gain,
        "out": args.out_tsv,
    })


if __name__ == "__main__":
    main()
