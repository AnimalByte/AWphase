#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if fields is None:
        fields = []
        seen = set()
        for r in rows:
            for k in r:
                if k not in seen:
                    fields.append(k)
                    seen.add(k)

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def state(row):
    for k in ("local_phase_state", "phase_state"):
        if k in row:
            try:
                return int(float(row.get(k) or 0))
            except Exception:
                return 0
    return 0

def usable_block(row):
    b = str(row.get("block_id", "")).strip()
    if b == "" or b.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return b

def load_calls(path):
    out = {}
    for r in read_tsv(path):
        try:
            p = int(float(r["pos"]))
        except Exception:
            continue
        out[p] = {
            "pos": p,
            "state": state(r),
            "block_id": usable_block(r),
            "row": r,
        }
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-tsv", required=True, help="AWPhase local_calls.tsv")
    ap.add_argument("--reference-tsv", required=True, help="WhatsHap local_calls.tsv")
    ap.add_argument("--out-block-tsv", required=True)
    ap.add_argument("--out-site-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    ap.add_argument("--min-overlap-sites", type=int, default=2)
    args = ap.parse_args()

    target = load_calls(args.target_tsv)
    ref = load_calls(args.reference_tsv)

    by_block = defaultdict(list)

    for pos, t in target.items():
        ts = t["state"]
        tb = t["block_id"]
        if ts == 0 or not tb:
            continue
        r = ref.get(pos)
        if not r:
            continue
        rs = r["state"]
        if rs == 0:
            continue
        by_block[tb].append((pos, ts, rs))

    block_rows = []
    site_rows = []

    total_overlap = 0
    total_residual_errors = 0
    total_pair_edges = 0
    total_pair_mismatches = 0
    mixed_blocks = 0
    pure_flipped_blocks = 0
    pure_kept_blocks = 0

    for block_id, items in sorted(by_block.items(), key=lambda kv: min(x[0] for x in kv[1])):
        items = sorted(items)
        same = sum(1 for _, ts, rs in items if ts == rs)
        opposite = sum(1 for _, ts, rs in items if ts == -rs)
        n = same + opposite

        if n < args.min_overlap_sites:
            continue

        # Best block orientation relative to reference.
        # +1 means keep AWPhase states; -1 means flip AWPhase states.
        orientation = 1 if same >= opposite else -1
        residual = min(same, opposite)
        residual_frac = residual / n if n else 0.0
        majority = max(same, opposite)

        if residual == 0 and orientation == -1:
            pure_flipped_blocks += 1
        elif residual == 0 and orientation == 1:
            pure_kept_blocks += 1
        elif residual > 0:
            mixed_blocks += 1

        # Adjacent relative-phase mismatch inside this AWPhase block.
        # This is invariant to whole-block flips.
        pair_edges = 0
        pair_mismatches = 0
        for (p1, ts1, rs1), (p2, ts2, rs2) in zip(items, items[1:]):
            aw_rel = ts1 * ts2
            ref_rel = rs1 * rs2
            pair_edges += 1
            if aw_rel != ref_rel:
                pair_mismatches += 1

        pair_mismatch_rate = pair_mismatches / pair_edges if pair_edges else 0.0

        total_overlap += n
        total_residual_errors += residual
        total_pair_edges += pair_edges
        total_pair_mismatches += pair_mismatches

        block_rows.append({
            "block_id": block_id,
            "start_pos": items[0][0],
            "end_pos": items[-1][0],
            "span_bp": items[-1][0] - items[0][0] + 1,
            "n_overlap": n,
            "same_count": same,
            "opposite_count": opposite,
            "best_orientation": orientation,
            "majority_agree": majority,
            "residual_errors_after_best_orientation": residual,
            "residual_error_frac": f"{residual_frac:.6f}",
            "adjacent_pair_edges": pair_edges,
            "adjacent_pair_mismatches": pair_mismatches,
            "adjacent_pair_mismatch_rate": f"{pair_mismatch_rate:.6f}",
        })

        for pos, ts, rs in items:
            oriented_ts = ts * orientation
            site_rows.append({
                "pos": pos,
                "block_id": block_id,
                "awphase_state": ts,
                "reference_state": rs,
                "best_block_orientation": orientation,
                "oriented_awphase_state": oriented_ts,
                "matches_after_best_block_orientation": int(oriented_ts == rs),
            })

    block_rows.sort(
        key=lambda r: (
            int(r["residual_errors_after_best_orientation"]),
            float(r["adjacent_pair_mismatch_rate"]),
            int(r["n_overlap"]),
        ),
        reverse=True,
    )

    summary = {
        "target_tsv": args.target_tsv,
        "reference_tsv": args.reference_tsv,
        "blocks_evaluated": len(block_rows),
        "total_overlap_sites": total_overlap,
        "total_residual_errors_after_best_block_orientation": total_residual_errors,
        "residual_error_rate_after_best_block_orientation": (
            total_residual_errors / total_overlap if total_overlap else None
        ),
        "total_adjacent_pair_edges": total_pair_edges,
        "total_adjacent_pair_mismatches": total_pair_mismatches,
        "adjacent_pair_mismatch_rate": (
            total_pair_mismatches / total_pair_edges if total_pair_edges else None
        ),
        "pure_kept_blocks": pure_kept_blocks,
        "pure_flipped_blocks": pure_flipped_blocks,
        "mixed_blocks": mixed_blocks,
    }

    write_tsv(args.out_block_tsv, block_rows)
    write_tsv(args.out_site_tsv, site_rows)
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
