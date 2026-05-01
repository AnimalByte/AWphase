#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path
from bisect import bisect_left
from collections import defaultdict

import numpy as np
import pysam

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

def safe_log1p(x):
    x = fval(x, 0.0)
    if x < 0:
        x = 0.0
    return math.log1p(x)

def entropy_binary(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))

def load_variants(path):
    obj = json.load(open(path))
    rows = obj.get("variants", obj) if isinstance(obj, dict) else obj

    out = {}
    for r in rows:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        ref = str(r.get("ref_allele", r.get("ref", ""))).upper()
        alt = str(r.get("alt_allele", r.get("alt", ""))).upper()
        out[pos] = (ref, alt)
    return out

def load_genetic_map(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    bp = []
    cm = []

    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                cm_val = float(parts[-2])
                bp_val = int(float(parts[-1]))
            except Exception:
                continue
            bp.append(bp_val)
            cm.append(cm_val)

    if not bp:
        raise SystemExit(f"No genetic-map points loaded from {path}")

    order = np.argsort(bp)
    bp = np.array([bp[i] for i in order], dtype=np.float64)
    cm = np.array([cm[i] for i in order], dtype=np.float64)
    return bp, cm

def interp_cm(pos, map_bp, map_cm):
    return float(np.interp(float(pos), map_bp, map_cm))

def load_phase6_calls(path):
    rows = read_tsv(path)

    anchors_by_block = defaultdict(list)
    rows_by_pos = {}

    for r in rows:
        pos = ival(r.get("pos"), None)
        if pos is None:
            continue

        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))
        rows_by_pos[pos] = r

        if state != 0 and block not in {"", "0", "unassigned"}:
            anchors_by_block[block].append((pos, state))

    for b in anchors_by_block:
        anchors_by_block[b].sort(key=lambda x: x[0])

    return rows, rows_by_pos, anchors_by_block

def load_truth_and_block_orientation(site_comparison_tsv):
    rows = read_tsv(site_comparison_tsv)

    truth_by_pos = {}
    block_counts = defaultdict(lambda: [0, 0])  # same, opposite

    for r in rows:
        pos = ival(r.get("pos"), None)
        if pos is None:
            continue

        pred = ival(r.get("pred_state", r.get("local_phase_state", 0)), 0)
        truth = ival(r.get("truth_state", 0), 0)
        block = str(r.get("block_id", "unassigned"))

        if truth != 0:
            truth_by_pos[pos] = truth

        if pred != 0 and truth != 0 and block not in {"", "0", "unassigned"}:
            if pred == truth:
                block_counts[block][0] += 1
            elif -pred == truth:
                block_counts[block][1] += 1

    orientation = {}
    orientation_margin = {}

    for block, (same, opp) in block_counts.items():
        if same > opp:
            orientation[block] = 1
            orientation_margin[block] = same - opp
        elif opp > same:
            orientation[block] = -1
            orientation_margin[block] = opp - same
        else:
            orientation[block] = 0
            orientation_margin[block] = 0

    return truth_by_pos, orientation, orientation_margin

def block_stats(anchors_by_block):
    out = {}
    for block, arr in anchors_by_block.items():
        poss = [x[0] for x in arr]
        if not poss:
            continue
        span = max(poss) - min(poss) + 1
        out[block] = {
            "scaffold_block_anchor_count": len(poss),
            "scaffold_block_span_bp": span,
            "scaffold_block_density_per_kb": len(poss) / max(span / 1000.0, 1e-9),
        }
    return out

def choose_scaffold_block(pos, anchors_by_block, map_bp, map_cm, max_cm, max_each_side, min_anchors):
    qcm = interp_cm(pos, map_bp, map_cm)

    best = None

    for block, arr in anchors_by_block.items():
        if not arr:
            continue

        poss = [x[0] for x in arr]
        cms = [interp_cm(p, map_bp, map_cm) for p in poss]
        i = bisect_left(poss, pos)

        chosen_left = []
        j = i - 1
        while j >= 0 and abs(qcm - cms[j]) <= max_cm and len(chosen_left) < max_each_side:
            chosen_left.append((arr[j][0], arr[j][1], abs(qcm - cms[j]), abs(pos - arr[j][0])))
            j -= 1

        chosen_right = []
        j = i
        while j < len(arr) and abs(cms[j] - qcm) <= max_cm and len(chosen_right) < max_each_side:
            chosen_right.append((arr[j][0], arr[j][1], abs(cms[j] - qcm), abs(arr[j][0] - pos)))
            j += 1

        chosen = list(reversed(chosen_left)) + chosen_right
        if len(chosen) < min_anchors:
            continue

        two_sided = 1 if chosen_left and chosen_right else 0
        mean_cm_dist = sum(x[2] for x in chosen) / len(chosen)
        min_cm_dist = min(x[2] for x in chosen)

        # Prefer two-sided scaffolds, then more anchors, then closer anchors.
        score = (
            two_sided,
            len(chosen),
            -mean_cm_dist,
            -min_cm_dist,
        )

        if best is None or score > best["score"]:
            best = {
                "block": block,
                "anchors": [(x[0], x[1]) for x in chosen],
                "left_count": len(chosen_left),
                "right_count": len(chosen_right),
                "two_sided": two_sided,
                "mean_cm_dist": mean_cm_dist,
                "min_cm_dist": min_cm_dist,
                "max_cm_dist": max(x[2] for x in chosen),
                "mean_bp_dist": sum(x[3] for x in chosen) / len(chosen),
                "min_bp_dist": min(x[3] for x in chosen),
                "max_bp_dist": max(x[3] for x in chosen),
                "score": score,
            }

    return best

def load_panel_haps(panel_bcf, chrom, start, end, needed_variants):
    vf = pysam.VariantFile(panel_bcf)
    samples = list(vf.header.samples)
    n_haps = len(samples) * 2

    states = {}
    found = 0
    allele_mismatch = 0
    nonbiallelic = 0

    needed = set(needed_variants)

    for rec in vf.fetch(chrom, max(0, start - 1), end):
        pos = int(rec.pos)
        if pos not in needed:
            continue

        want_ref, want_alt = needed_variants[pos]
        ref = str(rec.ref).upper()
        alts = rec.alts or []

        if len(alts) != 1:
            nonbiallelic += 1
            continue

        alt = str(alts[0]).upper()
        if ref != want_ref or alt != want_alt:
            allele_mismatch += 1
            continue

        alleles = np.zeros(n_haps, dtype=np.int8)
        valid = np.zeros(n_haps, dtype=bool)

        for si, sample in enumerate(samples):
            call = rec.samples[sample]
            gt = call.get("GT")
            phased = bool(call.phased)

            if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None:
                continue
            if not phased:
                continue
            if gt[0] not in (0, 1) or gt[1] not in (0, 1):
                continue

            alleles[2 * si] = gt[0]
            alleles[2 * si + 1] = gt[1]
            valid[2 * si] = True
            valid[2 * si + 1] = True

        states[pos] = (alleles, valid)
        found += 1

    vf.close()

    return {
        "samples": samples,
        "n_haps": n_haps,
        "states": states,
        "positions_needed": len(needed),
        "positions_found": found,
        "allele_mismatch": allele_mismatch,
        "nonbiallelic": nonbiallelic,
    }

def effective_n(weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w[w > 0]
    if w.size == 0:
        return 0.0
    return float((w.sum() ** 2) / max(float((w * w).sum()), 1e-12))

def topk_scores(score, valid, k):
    idx = np.where(valid)[0]
    if idx.size == 0:
        return idx, np.array([], dtype=np.float64)

    vals = score[idx]
    if vals.size > k:
        part = np.argpartition(vals, -k)[-k:]
        idx = idx[part]
        vals = vals[part]

    order = np.argsort(vals)[::-1]
    return idx[order], vals[order]

def projection_features(pos, anchors, panel, map_bp, map_cm, decay_cm, top_k):
    empty = {
        "phase8d_candidate_in_panel": 0,
        "phase8d_panel_usable_anchors": 0,
        "phase8d_support_plus": 0.0,
        "phase8d_support_minus": 0.0,
        "phase8d_total_support": 0.0,
        "phase8d_margin": 0.0,
        "phase8d_confidence": 0.0,
        "phase8d_entropy": 0.0,
        "phase8d_top_plus_score": 0.0,
        "phase8d_top_minus_score": 0.0,
        "phase8d_top_margin": 0.0,
        "phase8d_effective_donors": 0.0,
        "phase8d_hmm_log_margin": 0.0,
        "phase8d_abs_hmm_log_margin": 0.0,
        "pred_state": 0,
    }

    states = panel["states"]
    if pos not in states:
        return dict(empty)

    cand_alleles, cand_valid = states[pos]
    n_haps = panel["n_haps"]

    plus_score = np.zeros(n_haps, dtype=np.float64)
    minus_score = np.zeros(n_haps, dtype=np.float64)
    plus_valid = np.zeros(n_haps, dtype=bool)
    minus_valid = np.zeros(n_haps, dtype=bool)

    qcm = interp_cm(pos, map_bp, map_cm)
    usable = 0

    for apos, astate in anchors:
        if apos not in states:
            continue

        a_alleles, a_valid = states[apos]
        acm = interp_cm(apos, map_bp, map_cm)
        dcm = abs(qcm - acm)
        w = math.exp(-dcm / max(decay_cm, 1e-9))

        plus_expected = 1 if astate == 1 else 0
        minus_expected = 1 - plus_expected

        plus_match = a_valid & (a_alleles == plus_expected)
        minus_match = a_valid & (a_alleles == minus_expected)

        plus_score[plus_match] += w
        minus_score[minus_match] += w
        plus_valid |= plus_match
        minus_valid |= minus_match
        usable += 1

    if usable == 0:
        out = dict(empty)
        out["phase8d_candidate_in_panel"] = 1
        return out

    plus_idx, plus_vals = topk_scores(plus_score, plus_valid & cand_valid, top_k)
    minus_idx, minus_vals = topk_scores(minus_score, minus_valid & cand_valid, top_k)

    support_plus = 0.0
    support_minus = 0.0

    # Candidate state +1 means candidate goes on same local haplotype as scaffold plus.
    for idx, sc in zip(plus_idx, plus_vals):
        if cand_alleles[idx] == 1:
            support_plus += float(sc)
        else:
            support_minus += float(sc)

    for idx, sc in zip(minus_idx, minus_vals):
        if cand_alleles[idx] == 0:
            support_plus += float(sc)
        else:
            support_minus += float(sc)

    total = support_plus + support_minus
    if total > 0:
        major = max(support_plus, support_minus)
        minor = min(support_plus, support_minus)
        conf = major / total
        log_margin = math.log((support_plus + 1.0) / (support_minus + 1.0))
        pred_state = 1 if support_plus >= support_minus else -1
    else:
        major = minor = conf = log_margin = 0.0
        pred_state = 0

    all_top = list(plus_vals) + list(minus_vals)
    top_plus = float(plus_vals[0]) if len(plus_vals) else 0.0
    top_minus = float(minus_vals[0]) if len(minus_vals) else 0.0

    out = dict(empty)
    out["phase8d_candidate_in_panel"] = 1
    out["phase8d_panel_usable_anchors"] = usable
    out["phase8d_support_plus"] = support_plus
    out["phase8d_support_minus"] = support_minus
    out["phase8d_total_support"] = total
    out["phase8d_margin"] = major - minor
    out["phase8d_confidence"] = conf
    out["phase8d_entropy"] = entropy_binary(conf)
    out["phase8d_top_plus_score"] = top_plus
    out["phase8d_top_minus_score"] = top_minus
    out["phase8d_top_margin"] = abs(top_plus - top_minus)
    out["phase8d_effective_donors"] = effective_n(all_top)
    out["phase8d_hmm_log_margin"] = log_margin
    out["phase8d_abs_hmm_log_margin"] = abs(log_margin)
    out["log_phase8d_margin"] = safe_log1p(out["phase8d_margin"])
    out["log_phase8d_total_support"] = safe_log1p(total)
    out["log_phase8d_effective_donors"] = safe_log1p(out["phase8d_effective_donors"])
    out["pred_state"] = pred_state
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-window", required=True)
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--truth-comparison-phase6c-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--genetic-map", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--max-cm", type=float, default=0.25)
    ap.add_argument("--decay-cm", type=float, default=0.05)
    ap.add_argument("--max-anchors-each-side", type=int, default=16)
    ap.add_argument("--min-anchors", type=int, default=2)
    ap.add_argument("--top-k", type=int, default=64)
    ap.add_argument("--min-confidence-rule", type=float, default=0.85)
    ap.add_argument("--min-margin-rule", type=float, default=2.0)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    variants = load_variants(args.variant_json)
    local_rows, rows_by_pos, anchors_by_block = load_phase6_calls(args.local_calls_phase6c_tsv)
    truth_by_pos, block_orientation, block_orientation_margin = load_truth_and_block_orientation(args.truth_comparison_phase6c_tsv)
    stats = block_stats(anchors_by_block)
    map_bp, map_cm = load_genetic_map(args.genetic_map)

    candidates = []

    for r in local_rows:
        pos = ival(r.get("pos"), None)
        if pos is None:
            continue

        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        if state != 0:
            continue

        if pos not in variants:
            continue

        scaffold = choose_scaffold_block(
            pos,
            anchors_by_block,
            map_bp,
            map_cm,
            args.max_cm,
            args.max_anchors_each_side,
            args.min_anchors,
        )

        if scaffold is None:
            continue

        candidates.append((pos, scaffold))

    needed_positions = set()
    for pos, scaffold in candidates:
        needed_positions.add(pos)
        for apos, _ in scaffold["anchors"]:
            needed_positions.add(apos)

    needed_variants = {p: variants[p] for p in needed_positions if p in variants}

    panel = load_panel_haps(
        args.panel_bcf,
        args.chrom,
        args.start,
        args.end,
        needed_variants,
    )

    out_rows = []

    for pos, scaffold in candidates:
        block = scaffold["block"]
        anchors = scaffold["anchors"]

        feats = projection_features(
            pos,
            anchors,
            panel,
            map_bp,
            map_cm,
            args.decay_cm,
            args.top_k,
        )

        pred_state = ival(feats.get("pred_state"), 0)
        truth_state = truth_by_pos.get(pos, 0)
        orient = block_orientation.get(block, 0)

        if pred_state == 0:
            continue
        if truth_state == 0:
            continue
        if orient == 0:
            continue

        label = 1 if pred_state * orient == truth_state else 0

        bs = stats.get(block, {})

        row = {
            "source_window": args.source_window,
            "pos": pos,
            "selected_block_id": block,
            "pred_state": pred_state,
            "truth_state": truth_state,
            "label": label,

            "scaffold_anchor_count": len(anchors),
            "scaffold_left_anchor_count": scaffold["left_count"],
            "scaffold_right_anchor_count": scaffold["right_count"],
            "scaffold_two_sided": scaffold["two_sided"],
            "scaffold_mean_cm_dist": scaffold["mean_cm_dist"],
            "scaffold_min_cm_dist": scaffold["min_cm_dist"],
            "scaffold_max_cm_dist": scaffold["max_cm_dist"],
            "scaffold_mean_bp_dist": scaffold["mean_bp_dist"],
            "scaffold_min_bp_dist": scaffold["min_bp_dist"],
            "scaffold_max_bp_dist": scaffold["max_bp_dist"],

            "scaffold_block_anchor_count": bs.get("scaffold_block_anchor_count", 0),
            "scaffold_block_span_bp": bs.get("scaffold_block_span_bp", 0),
            "scaffold_block_density_per_kb": bs.get("scaffold_block_density_per_kb", 0.0),
            "log_scaffold_block_anchor_count": safe_log1p(bs.get("scaffold_block_anchor_count", 0)),
            "log_scaffold_block_span_bp": safe_log1p(bs.get("scaffold_block_span_bp", 0)),

            "phase8d_rule_accept": 1 if (
                feats.get("phase8d_candidate_in_panel", 0) == 1
                and feats.get("phase8d_confidence", 0.0) >= args.min_confidence_rule
                and feats.get("phase8d_margin", 0.0) >= args.min_margin_rule
                and len(anchors) >= args.min_anchors
            ) else 0,
        }

        row.update(feats)
        out_rows.append(row)

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

    positives = sum(int(r["label"]) for r in out_rows)
    negatives = len(out_rows) - positives

    print(json.dumps({
        "source_window": args.source_window,
        "phase6_unphased_candidates_considered": len(candidates),
        "training_rows": len(out_rows),
        "positives": positives,
        "negatives": negatives,
        "panel_positions_needed": panel["positions_needed"],
        "panel_positions_found": panel["positions_found"],
        "panel_allele_mismatch": panel["allele_mismatch"],
        "panel_nonbiallelic": panel["nonbiallelic"],
        "out": args.out_tsv,
    }, indent=2))

if __name__ == "__main__":
    main()
