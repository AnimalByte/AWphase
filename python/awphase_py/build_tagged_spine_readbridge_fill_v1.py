#!/usr/bin/env python3
import argparse
import bisect
import csv
import json
from collections import defaultdict
from pathlib import Path

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

def state_of(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def set_state(row, state):
    row["local_phase_state"] = int(state)
    row["phase_state"] = int(state)

def weight_of_obs(row):
    for k in ("weight_v3", "weight_v2", "allele_confidence", "effective_baseq", "raw_baseq"):
        if k in row and str(row.get(k, "")).strip() != "":
            v = fval(row.get(k), 0.0)
            if k in ("effective_baseq", "raw_baseq"):
                return max(0.1, min(3.0, v / 15.0))
            return max(0.1, min(6.0, v))
    return 1.0

def interval_id_for_pos(pos, anchors):
    i = bisect.bisect_left(anchors, pos)
    left = anchors[i - 1] if i > 0 else None
    right = anchors[i] if i < len(anchors) else None
    return (left, right)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tagged-local-calls-tsv", required=True)
    ap.add_argument("--fill-local-calls-tsv", required=True)
    ap.add_argument("--obs-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)

    ap.add_argument("--max-anchor-dist-bp", type=int, default=75000)
    ap.add_argument("--min-bridge-pairs", type=int, default=2)
    ap.add_argument("--min-orientation-margin", type=float, default=2.0)
    ap.add_argument("--min-obs-weight", type=float, default=0.2)
    ap.add_argument("--allow-one-sided-anchors", action="store_true")
    ap.add_argument("--fill-unbridged-as-missing", action="store_true")
    args = ap.parse_args()

    tagged_rows = read_tsv(args.tagged_local_calls_tsv)
    fill_rows = read_tsv(args.fill_local_calls_tsv)

    tagged_by_pos = {ival(r["pos"]): r for r in tagged_rows}
    fill_by_pos = {ival(r["pos"]): r for r in fill_rows}

    anchor_state = {}
    for pos, r in tagged_by_pos.items():
        st = state_of(r)
        if st != 0:
            anchor_state[pos] = st

    anchors = sorted(anchor_state)

    # Fill candidates are positions where tagged abstains but fill source has a nonzero call.
    fill_state = {}
    for pos, fr in fill_by_pos.items():
        if pos not in tagged_by_pos:
            continue
        if state_of(tagged_by_pos[pos]) != 0:
            continue
        st = state_of(fr)
        if st != 0:
            fill_state[pos] = st

    fill_positions = set(fill_state)
    anchor_positions = set(anchor_state)

    interval_for_fill = {}
    interval_fill_positions = defaultdict(list)
    for pos in sorted(fill_positions):
        iid = interval_id_for_pos(pos, anchors)
        left, right = iid

        if not args.allow_one_sided_anchors and (left is None or right is None):
            continue

        # Require at least one nearby anchor, or both if available.
        near_left = left is not None and abs(pos - left) <= args.max_anchor_dist_bp
        near_right = right is not None and abs(right - pos) <= args.max_anchor_dist_bp
        if not (near_left or near_right):
            continue

        interval_for_fill[pos] = iid
        interval_fill_positions[iid].append(pos)

    # Read observations for anchor and fill positions.
    # Encoding consistency rule:
    # For positions i,j on the same read, predicted phase states s_i,s_j and observed alleles a_i,a_j
    # are consistent if a_i*s_i == a_j*s_j.
    obs_by_read = defaultdict(dict)

    with open(args.obs_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rid = str(r.get("read_id", "")).strip()
            if not rid:
                continue
            pos = ival(r.get("site_pos", r.get("pos", "")), None)
            if pos is None:
                continue
            if pos not in anchor_positions and pos not in fill_positions:
                continue

            allele = ival(r.get("allele", 0), 0)
            if allele not in (-1, 1):
                continue
            if ival(r.get("is_ambiguous", 0), 0) != 0:
                continue

            w = weight_of_obs(r)
            if w < args.min_obs_weight:
                continue

            obs_by_read[rid][pos] = (allele, w)

    interval_score = defaultdict(lambda: {
        "keep": 0.0,
        "flip": 0.0,
        "pairs": 0,
        "reads": set(),
        "anchors_seen": set(),
        "fills_seen": set(),
    })

    for rid, obs in obs_by_read.items():
        read_anchor_positions = [p for p in obs if p in anchor_positions]
        read_fill_positions = [p for p in obs if p in fill_positions and p in interval_for_fill]

        if not read_anchor_positions or not read_fill_positions:
            continue

        for fp in read_fill_positions:
            iid = interval_for_fill[fp]
            f_allele, f_w = obs[fp]
            f_state = fill_state[fp]

            for ap in read_anchor_positions:
                if abs(fp - ap) > args.max_anchor_dist_bp:
                    continue

                # Only use anchors belonging to the same bracketing interval if possible.
                left, right = iid
                if left is not None and right is not None and ap not in (left, right):
                    continue

                a_allele, a_w = obs[ap]
                a_state = anchor_state[ap]
                w = min(f_w, a_w)

                keep_consistent = (a_allele * a_state) == (f_allele * f_state)
                flip_consistent = (a_allele * a_state) == (f_allele * (-f_state))

                interval_score[iid]["keep"] += w if keep_consistent else -w
                interval_score[iid]["flip"] += w if flip_consistent else -w
                interval_score[iid]["pairs"] += 1
                interval_score[iid]["reads"].add(rid)
                interval_score[iid]["anchors_seen"].add(ap)
                interval_score[iid]["fills_seen"].add(fp)

    interval_orientation = {}
    interval_rows = []

    for iid, positions in sorted(interval_fill_positions.items(), key=lambda kv: (kv[0][0] or -1, kv[0][1] or 10**18)):
        sc = interval_score.get(iid, None)
        left, right = iid

        if sc is None:
            keep = flip = 0.0
            pairs = 0
            n_reads = 0
            n_anchors = 0
            n_fills = 0
        else:
            keep = sc["keep"]
            flip = sc["flip"]
            pairs = sc["pairs"]
            n_reads = len(sc["reads"])
            n_anchors = len(sc["anchors_seen"])
            n_fills = len(sc["fills_seen"])

        margin = abs(keep - flip)
        orientation = 0
        accepted = False
        reason = "no_bridge"

        if pairs >= args.min_bridge_pairs and margin >= args.min_orientation_margin:
            orientation = 1 if keep >= flip else -1
            accepted = True
            reason = "accepted_keep" if orientation == 1 else "accepted_flip"
        elif pairs < args.min_bridge_pairs:
            reason = "too_few_bridge_pairs"
        else:
            reason = "low_orientation_margin"

        interval_orientation[iid] = orientation if accepted else 0

        interval_rows.append({
            "left_anchor_pos": "" if left is None else left,
            "right_anchor_pos": "" if right is None else right,
            "left_anchor_state": 0 if left is None else anchor_state.get(left, 0),
            "right_anchor_state": 0 if right is None else anchor_state.get(right, 0),
            "n_fill_positions": len(positions),
            "keep_score": f"{keep:.6f}",
            "flip_score": f"{flip:.6f}",
            "orientation_margin": f"{margin:.6f}",
            "bridge_pairs": pairs,
            "bridge_reads": n_reads,
            "anchors_seen": n_anchors,
            "fills_seen": n_fills,
            "orientation": orientation if accepted else 0,
            "accepted": int(accepted),
            "reason": reason,
            "positions": ",".join(map(str, positions)),
        })

    # Build output from tagged backbone.
    out = [dict(r) for r in tagged_rows]
    by_pos_out = {ival(r["pos"]): r for r in out}

    filled = 0
    changed = 0
    skipped_unoriented = 0
    skipped_no_interval = 0

    for pos, st in fill_state.items():
        row = by_pos_out.get(pos)
        if row is None:
            continue

        old = state_of(row)
        if old != 0:
            continue

        iid = interval_for_fill.get(pos)
        if iid is None:
            skipped_no_interval += 1
            continue

        orient = interval_orientation.get(iid, 0)
        if orient == 0:
            skipped_unoriented += 1
            if args.fill_unbridged_as_missing:
                continue
            else:
                continue

        new_state = st * orient
        set_state(row, new_state)

        row["phase5_source"] = "tagged_spine_readbridge_fill"
        row["phase5_original_fill_state"] = st
        row["phase5_interval_orientation"] = orient
        row["phase5_interval_left"] = "" if iid[0] is None else iid[0]
        row["phase5_interval_right"] = "" if iid[1] is None else iid[1]

        filled += 1
        if old != new_state:
            changed += 1

    for r in out:
        r.setdefault("phase5_source", "tagged_backbone" if state_of(r) != 0 else "missing")

    write_tsv(args.out_tsv, out)

    summary = {
        "tagged_backbone_sites": len(anchor_state),
        "fill_candidates": len(fill_state),
        "intervals": len(interval_rows),
        "accepted_intervals": sum(ival(r["accepted"]) for r in interval_rows),
        "filled_sites": filled,
        "changed_sites": changed,
        "skipped_unoriented": skipped_unoriented,
        "skipped_no_interval": skipped_no_interval,
        "max_anchor_dist_bp": args.max_anchor_dist_bp,
        "min_bridge_pairs": args.min_bridge_pairs,
        "min_orientation_margin": args.min_orientation_margin,
    }

    out_summary = Path(args.out_summary_json)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")

    interval_path = out_summary.with_suffix(".intervals.tsv")
    write_tsv(interval_path, interval_rows)

    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
