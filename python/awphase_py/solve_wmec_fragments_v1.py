#!/usr/bin/env python3
import argparse
import csv
import json
import heapq
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def read_variant_json(path):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj

    rows = []
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        rows.append({
            "pos": pos,
            "ref": str(r.get("ref_allele", r.get("ref", ""))),
            "alt": str(r.get("alt_allele", r.get("alt", ""))),
        })

    rows.sort(key=lambda x: x["pos"])
    return rows

def load_fragments(path):
    frags = []
    for r in read_tsv(path):
        pos = json.loads(r["positions_json"])
        alleles = json.loads(r["alleles_json"])
        weights = json.loads(r["weights_json"])

        if len(pos) != len(alleles) or len(pos) != len(weights):
            continue

        items = []
        for p, a, w in zip(pos, alleles, weights):
            try:
                p = int(p)
                a = int(a)
                w = float(w)
            except Exception:
                continue
            if a not in {-1, 1}:
                continue
            if w <= 0:
                continue
            items.append((p, a, w))

        if len(items) >= 2:
            items.sort()
            frags.append({
                "fragment_id": r.get("fragment_id", ""),
                "items": items,
                "score": sum(w for _, _, w in items),
            })

    return frags

class DSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

def fragment_cost_local(x, sites, alleles, weights):
    plus = 0.0
    minus = 0.0
    for idx, a, w in zip(sites, alleles, weights):
        if x[idx] != a:
            plus += w
        if x[idx] != -a:
            minus += w
    return min(plus, minus)

def fragment_orientation(x, sites, alleles, weights):
    plus = 0.0
    minus = 0.0
    for idx, a, w in zip(sites, alleles, weights):
        if x[idx] != a:
            plus += w
        if x[idx] != -a:
            minus += w
    return 1 if plus <= minus else -1

def total_cost(x, frags):
    return sum(fragment_cost_local(x, f["sites"], f["alleles"], f["weights"]) for f in frags)

def exact_solve(n, frags):
    # Fix first site to +1 because global polarity is arbitrary.
    best = None
    second = None
    best_x = None

    if n == 0:
        return [], 0.0, 0.0

    for mask in range(1 << max(0, n - 1)):
        x = [1] * n
        for i in range(1, n):
            x[i] = 1 if ((mask >> (i - 1)) & 1) else -1

        c = total_cost(x, frags)

        if best is None or c < best:
            second = best
            best = c
            best_x = x
        elif second is None or c < second:
            second = c

    margin = 0.0 if second is None else float(second - best)
    return best_x, float(best), margin

def greedy_init(n, frags):
    pair = defaultdict(float)

    for f in frags:
        sites = f["sites"]
        alleles = f["alleles"]
        weights = f["weights"]
        k = len(sites)

        # Reads usually cover few hets. Cap all-pairs for safety.
        if k > 20:
            continue

        for i in range(k):
            for j in range(i + 1, k):
                a = sites[i]
                b = sites[j]
                sign = alleles[i] * alleles[j]
                w = min(weights[i], weights[j])
                if a > b:
                    a, b = b, a
                pair[(a, b)] += sign * w

    adj = defaultdict(list)
    for (a, b), v in pair.items():
        if v == 0:
            continue
        sign = 1 if v > 0 else -1
        w = abs(v)
        adj[a].append((b, sign, w))
        adj[b].append((a, sign, w))

    x = [0] * n

    for start in range(n):
        if x[start] != 0:
            continue

        x[start] = 1
        heap = []
        for nb, sign, w in adj[start]:
            heapq.heappush(heap, (-w, start, nb, sign))

        while heap:
            negw, src, dst, sign = heapq.heappop(heap)
            if x[dst] != 0:
                continue
            x[dst] = x[src] * sign
            for nb, sign2, w2 in adj[dst]:
                if x[nb] == 0:
                    heapq.heappush(heap, (-w2, dst, nb, sign2))

    for i in range(n):
        if x[i] == 0:
            x[i] = 1

    return x

def local_refine(x, frags, max_iter=5):
    by_site = defaultdict(list)
    costs = []

    for fi, f in enumerate(frags):
        c = fragment_cost_local(x, f["sites"], f["alleles"], f["weights"])
        costs.append(c)
        for s in f["sites"]:
            by_site[s].append(fi)

    best_cost = sum(costs)

    for _ in range(max_iter):
        improved = False

        for s in range(len(x)):
            affected = by_site.get(s, [])
            if not affected:
                continue

            old_val = x[s]
            old_sum = sum(costs[fi] for fi in affected)

            x[s] = -x[s]
            new_costs = {}
            new_sum = 0.0

            for fi in affected:
                f = frags[fi]
                c = fragment_cost_local(x, f["sites"], f["alleles"], f["weights"])
                new_costs[fi] = c
                new_sum += c

            delta = new_sum - old_sum

            if delta < -1e-9:
                for fi, c in new_costs.items():
                    costs[fi] = c
                best_cost += delta
                improved = True
            else:
                x[s] = old_val

        if not improved:
            break

    return x, float(best_cost)

def solve_greedy(n, frags, max_iter):
    x = greedy_init(n, frags)
    x, c = local_refine(x, frags, max_iter=max_iter)
    return x, c, ""

def compute_site_support(x, frags):
    support = defaultdict(float)
    conflict = defaultdict(float)

    for f in frags:
        o = fragment_orientation(x, f["sites"], f["alleles"], f["weights"])
        for idx, a, w in zip(f["sites"], f["alleles"], f["weights"]):
            if o * a == x[idx]:
                support[idx] += w
            else:
                conflict[idx] += w

    return support, conflict

def project_fragments_to_positions(frags, positions):
    pos_set = set(positions)
    pos_to_idx = {p: i for i, p in enumerate(positions)}

    out = []
    for f in frags:
        sites = []
        alleles = []
        weights = []

        for p, a, w in f["items"]:
            if p in pos_set:
                sites.append(pos_to_idx[p])
                alleles.append(a)
                weights.append(w)

        if len(sites) >= 2:
            out.append({
                "fragment_id": f["fragment_id"],
                "sites": sites,
                "alleles": alleles,
                "weights": weights,
            })

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fragments-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--out-local-calls-tsv", required=True)
    ap.add_argument("--out-components-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--max-exact-sites", type=int, default=18)
    ap.add_argument("--max-component-sites", type=int, default=256)
    ap.add_argument("--min-component-sites", type=int, default=2)
    ap.add_argument("--min-fragments-per-component", type=int, default=1)
    ap.add_argument("--local-refine-iters", type=int, default=5)
    args = ap.parse_args()

    variant_rows = read_variant_json(args.variant_json)
    variant_positions = [r["pos"] for r in variant_rows]
    variant_pos_set = set(variant_positions)

    raw_frags = load_fragments(args.fragments_tsv)
    raw_frags = [
        {
            "fragment_id": f["fragment_id"],
            "items": [(p, a, w) for p, a, w in f["items"] if p in variant_pos_set],
            "score": f["score"],
        }
        for f in raw_frags
    ]
    raw_frags = [f for f in raw_frags if len(f["items"]) >= 2]

    dsu = DSU()
    for f in raw_frags:
        ps = [p for p, _, _ in f["items"]]
        for p in ps:
            dsu.find(p)
        for p in ps[1:]:
            dsu.union(ps[0], p)

    comp_pos = defaultdict(set)
    for f in raw_frags:
        for p, _, _ in f["items"]:
            comp_pos[dsu.find(p)].add(p)

    raw_components = [sorted(v) for v in comp_pos.values() if len(v) >= args.min_component_sites]
    raw_components.sort(key=lambda x: (x[0], x[-1]))

    # Split large connected components into position chunks for a safe prototype.
    components = []
    for comp in raw_components:
        if len(comp) <= args.max_component_sites:
            components.append(comp)
        else:
            for i in range(0, len(comp), args.max_component_sites):
                chunk = comp[i:i + args.max_component_sites]
                if len(chunk) >= args.min_component_sites:
                    components.append(chunk)

    solved = {}
    component_rows = []

    exact_components = 0
    greedy_components = 0
    skipped_components = 0
    total_projected_frags = 0

    for ci, positions in enumerate(components, start=1):
        projected = project_fragments_to_positions(raw_frags, positions)

        if len(projected) < args.min_fragments_per_component:
            skipped_components += 1
            continue

        n = len(positions)
        total_projected_frags += len(projected)

        if n <= args.max_exact_sites:
            x, obj, margin = exact_solve(n, projected)
            solver = "exact"
            exact_components += 1
        else:
            x, obj, margin = solve_greedy(n, projected, args.local_refine_iters)
            solver = "greedy_refine"
            greedy_components += 1

        support, conflict = compute_site_support(x, projected)

        block_id = f"phase6_comp_{ci}"

        for idx, pos in enumerate(positions):
            sup = support.get(idx, 0.0)
            con = conflict.get(idx, 0.0)
            conf = sup / (sup + con) if (sup + con) > 0 else 0.0

            solved[pos] = {
                "pos": pos,
                "block_id": block_id,
                "phase_state": x[idx],
                "local_phase_state": x[idx],
                "confidence": f"{conf:.6f}",
                "phase6_solver": solver,
                "phase6_component_id": block_id,
                "phase6_component_sites": n,
                "phase6_component_fragments": len(projected),
                "phase6_objective": f"{obj:.6f}",
                "phase6_margin": margin,
                "phase6_site_support": f"{sup:.6f}",
                "phase6_site_conflict": f"{con:.6f}",
            }

        component_rows.append({
            "component_id": block_id,
            "start_pos": positions[0],
            "end_pos": positions[-1],
            "span_bp": positions[-1] - positions[0] + 1,
            "n_sites": n,
            "n_fragments": len(projected),
            "solver": solver,
            "objective": f"{obj:.6f}",
            "margin": margin,
        })

    out_rows = []
    for r in variant_rows:
        pos = r["pos"]
        if pos in solved:
            out_rows.append(solved[pos])
        else:
            out_rows.append({
                "pos": pos,
                "block_id": "unassigned",
                "phase_state": 0,
                "local_phase_state": 0,
                "confidence": "0.000000",
                "phase6_solver": "unassigned",
                "phase6_component_id": "",
                "phase6_component_sites": 0,
                "phase6_component_fragments": 0,
                "phase6_objective": "",
                "phase6_margin": "",
                "phase6_site_support": "0.000000",
                "phase6_site_conflict": "0.000000",
            })

    Path(args.out_local_calls_tsv).parent.mkdir(parents=True, exist_ok=True)

    local_fields = [
        "pos", "block_id", "phase_state", "local_phase_state", "confidence",
        "phase6_solver", "phase6_component_id", "phase6_component_sites",
        "phase6_component_fragments", "phase6_objective", "phase6_margin",
        "phase6_site_support", "phase6_site_conflict",
    ]

    with open(args.out_local_calls_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=local_fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    comp_fields = [
        "component_id", "start_pos", "end_pos", "span_bp", "n_sites",
        "n_fragments", "solver", "objective", "margin",
    ]

    with open(args.out_components_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=comp_fields, delimiter="\t")
        w.writeheader()
        w.writerows(component_rows)

    summary = {
        "fragments_tsv": args.fragments_tsv,
        "variant_json": args.variant_json,
        "raw_fragments": len(raw_frags),
        "raw_connected_components": len(raw_components),
        "components_after_chunking": len(components),
        "components_solved": len(component_rows),
        "components_skipped": skipped_components,
        "exact_components": exact_components,
        "greedy_components": greedy_components,
        "total_projected_fragments": total_projected_frags,
        "variant_rows": len(variant_rows),
        "solved_sites": len(solved),
        "out_local_calls_tsv": args.out_local_calls_tsv,
        "max_exact_sites": args.max_exact_sites,
        "max_component_sites": args.max_component_sites,
    }

    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
