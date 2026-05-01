#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict, deque
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


def get_relation(row):
    if "pred_relation" in row and str(row.get("pred_relation", "")).strip() != "":
        return ival(row.get("pred_relation"), 0)
    return ival(row.get("bridge_pred_relation"), 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--bridge-predictions-tsv", required=True)
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    rows = read_tsv(args.local_calls_tsv)
    preds = read_tsv(args.bridge_predictions_tsv)

    graph = defaultdict(list)
    accepted_edges = 0

    for r in preds:
        if str(r.get("source_window")) != args.source_window:
            continue

        p = fval(r.get("xgb_oof_probability"), -1.0)
        if p < args.threshold:
            continue

        a = str(r.get("block_a", ""))
        b = str(r.get("block_b", ""))
        rel = get_relation(r)

        if not a or not b or a == b or rel not in {-1, 1}:
            continue

        graph[a].append((b, rel))
        graph[b].append((a, rel))
        accepted_edges += 1

    orient = {}
    root_of = {}
    bad_roots = set()
    components = defaultdict(list)

    for start in sorted(graph):
        if start in orient:
            continue

        root = start
        orient[start] = 1
        root_of[start] = root
        components[root].append(start)
        q = deque([start])
        conflict = False

        while q:
            x = q.popleft()
            for y, rel in graph[x]:
                expected = orient[x] * rel

                if y not in orient:
                    orient[y] = expected
                    root_of[y] = root
                    components[root].append(y)
                    q.append(y)
                elif orient[y] != expected:
                    conflict = True

        if conflict:
            bad_roots.add(root)

    out_rows = []
    merged_sites = 0
    skipped_conflict_sites = 0

    for r in rows:
        rr = dict(r)
        block = str(rr.get("block_id", "unassigned"))
        state = ival(rr.get("local_phase_state", rr.get("phase_state", 0)), 0)

        rr["phase9a_bridge_merged"] = 0
        rr["phase9a_bridge_component_conflict_skipped"] = 0
        rr["phase9a_original_block_id"] = block

        if block in root_of and state != 0:
            root = root_of[block]

            if root in bad_roots:
                rr["phase9a_bridge_component_conflict_skipped"] = 1
                skipped_conflict_sites += 1
            else:
                mult = orient[block]
                new_state = state * mult

                if "local_phase_state" in rr:
                    rr["local_phase_state"] = new_state
                if "phase_state" in rr:
                    rr["phase_state"] = new_state

                rr["block_id"] = "phase9a_bridge_" + root
                rr["phase9a_bridge_merged"] = 1
                merged_sites += 1

        out_rows.append(rr)

    fields = list(rows[0].keys()) if rows else []
    for f in [
        "phase9a_bridge_merged",
        "phase9a_bridge_component_conflict_skipped",
        "phase9a_original_block_id",
    ]:
        if f not in fields:
            fields.append(f)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    print({
        "source_window": args.source_window,
        "threshold": args.threshold,
        "accepted_edges": accepted_edges,
        "graph_blocks": len(graph),
        "components": len(components),
        "conflicting_components_skipped": len(bad_roots),
        "merged_phased_sites": merged_sites,
        "skipped_conflict_sites": skipped_conflict_sites,
        "out": args.out_tsv,
    })


if __name__ == "__main__":
    main()
