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

def effective_n(weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w[w > 0]
    if w.size == 0:
        return 0.0
    return float((w.sum() ** 2) / max(float((w * w).sum()), 1e-12))

def softmax(vals):
    vals = np.asarray(vals, dtype=np.float64)
    if vals.size == 0:
        return vals
    m = np.max(vals)
    ex = np.exp(vals - m)
    s = np.sum(ex)
    if s <= 0:
        return np.ones_like(vals) / len(vals)
    return ex / s

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

def load_phase6_anchors(local_calls_tsv):
    rows = read_tsv(local_calls_tsv)
    anchors_by_block = defaultdict(list)

    for r in rows:
        pos = ival(r.get("pos"), None)
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))

        if pos is None:
            continue
        if state != 0 and block not in {"", "0", "unassigned"}:
            anchors_by_block[block].append((pos, state))

    for b in anchors_by_block:
        anchors_by_block[b].sort(key=lambda x: x[0])

    return anchors_by_block

def choose_anchors(pos, block, anchors_by_block, map_bp, map_cm, max_cm, max_each_side):
    arr = anchors_by_block.get(block, [])
    if not arr:
        return []

    poss = [x[0] for x in arr]
    cms = [interp_cm(p, map_bp, map_cm) for p in poss]
    qcm = interp_cm(pos, map_bp, map_cm)

    i = bisect_left(poss, pos)

    left = []
    j = i - 1
    while j >= 0 and abs(qcm - cms[j]) <= max_cm and len(left) < max_each_side:
        left.append((arr[j][0], arr[j][1], abs(qcm - cms[j]), abs(pos - arr[j][0])))
        j -= 1

    right = []
    j = i
    while j < len(arr) and abs(cms[j] - qcm) <= max_cm and len(right) < max_each_side:
        right.append((arr[j][0], arr[j][1], abs(cms[j] - qcm), abs(arr[j][0] - pos)))
        j += 1

    return list(reversed(left)) + right

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

def score_context(pos, anchors, panel, map_bp, map_cm, decay_cm, top_k):
    states = panel["states"]
    n_haps = panel["n_haps"]

    empty = {
        "usable_anchors": 0,
        "plus_prob": 0.0,
        "minus_prob": 0.0,
        "confidence": 0.0,
        "margin": 0.0,
        "entropy": 0.0,
        "top_plus_score": 0.0,
        "top_minus_score": 0.0,
        "top_score_margin": 0.0,
        "effective_haps": 0.0,
        "mean_cm_dist": -1.0,
        "min_cm_dist": -1.0,
        "max_cm_dist": -1.0,
        "mean_bp_dist": -1.0,
        "min_bp_dist": -1.0,
        "max_bp_dist": -1.0,
    }

    if pos not in states or not anchors:
        return empty

    cand_alleles, cand_valid = states[pos]

    plus_score = np.zeros(n_haps, dtype=np.float64)
    minus_score = np.zeros(n_haps, dtype=np.float64)
    plus_valid = np.zeros(n_haps, dtype=bool)
    minus_valid = np.zeros(n_haps, dtype=bool)

    qcm = interp_cm(pos, map_bp, map_cm)

    cm_dists = []
    bp_dists = []
    usable = 0

    for apos, astate, dcm_pre, dbp_pre in anchors:
        if apos not in states:
            continue

        a_alleles, a_valid = states[apos]
        acm = interp_cm(apos, map_bp, map_cm)
        dcm = abs(qcm - acm)
        dbp = abs(pos - apos)

        w = math.exp(-dcm / max(decay_cm, 1e-9))

        # Local phase convention:
        # state +1 means local plus hap carries ALT at anchor.
        # state -1 means local plus hap carries REF at anchor.
        plus_expected = 1 if astate == 1 else 0
        minus_expected = 1 - plus_expected

        plus_match = a_valid & (a_alleles == plus_expected)
        minus_match = a_valid & (a_alleles == minus_expected)

        plus_score[plus_match] += w
        minus_score[minus_match] += w

        plus_valid |= plus_match
        minus_valid |= minus_match

        usable += 1
        cm_dists.append(dcm)
        bp_dists.append(dbp)

    if usable == 0:
        return empty

    def topk(score, valid):
        idx = np.where(valid & cand_valid)[0]
        if idx.size == 0:
            return idx, np.array([], dtype=np.float64)
        vals = score[idx]
        if vals.size > top_k:
            part = np.argpartition(vals, -top_k)[-top_k:]
            idx = idx[part]
            vals = vals[part]
        order = np.argsort(vals)[::-1]
        return idx[order], vals[order]

    plus_idx, plus_vals = topk(plus_score, plus_valid)
    minus_idx, minus_vals = topk(minus_score, minus_valid)

    support_plus = 0.0
    support_minus = 0.0

    if plus_vals.size:
        plus_post = softmax(plus_vals)
        for idx, w in zip(plus_idx, plus_post):
            if cand_alleles[idx] == 1:
                support_plus += float(w)
            else:
                support_minus += float(w)

    if minus_vals.size:
        minus_post = softmax(minus_vals)
        for idx, w in zip(minus_idx, minus_post):
            if cand_alleles[idx] == 0:
                support_plus += float(w)
            else:
                support_minus += float(w)

    total = support_plus + support_minus

    if total > 0:
        p_plus = support_plus / total
        p_minus = support_minus / total
        conf = max(p_plus, p_minus)
        margin = abs(p_plus - p_minus)
    else:
        p_plus = p_minus = conf = margin = 0.0

    all_scores = list(plus_vals) + list(minus_vals)

    out = dict(empty)
    out["usable_anchors"] = usable
    out["plus_prob"] = p_plus
    out["minus_prob"] = p_minus
    out["confidence"] = conf
    out["margin"] = margin
    out["entropy"] = entropy_binary(conf)
    out["top_plus_score"] = float(plus_vals[0]) if plus_vals.size else 0.0
    out["top_minus_score"] = float(minus_vals[0]) if minus_vals.size else 0.0
    out["top_score_margin"] = abs(out["top_plus_score"] - out["top_minus_score"])
    out["effective_haps"] = effective_n(all_scores)
    out["mean_cm_dist"] = float(sum(cm_dists) / len(cm_dists))
    out["min_cm_dist"] = float(min(cm_dists))
    out["max_cm_dist"] = float(max(cm_dists))
    out["mean_bp_dist"] = float(sum(bp_dists) / len(bp_dists))
    out["min_bp_dist"] = float(min(bp_dists))
    out["max_bp_dist"] = float(max(bp_dists))

    return out

def make_features_for_candidate(pos, block, anchors, panel, map_bp, map_cm, decay_cm, top_k):
    states = panel["states"]

    empty = {
        "phase8c_candidate_in_panel": 0,
        "phase8c_anchor_count": len(anchors),
        "phase8c_left_anchor_count": 0,
        "phase8c_right_anchor_count": 0,
        "phase8c_two_sided": 0,
        "phase8c_copy_plus_prob": 0.0,
        "phase8c_copy_minus_prob": 0.0,
        "phase8c_copy_confidence": 0.0,
        "phase8c_copy_margin": 0.0,
        "phase8c_copy_entropy": 0.0,
        "phase8c_copy_effective_haps": 0.0,
        "phase8c_top_plus_score": 0.0,
        "phase8c_top_minus_score": 0.0,
        "phase8c_top_score_margin": 0.0,
        "phase8c_left_plus_prob": 0.0,
        "phase8c_right_plus_prob": 0.0,
        "phase8c_left_right_agree": 0,
        "phase8c_left_confidence": 0.0,
        "phase8c_right_confidence": 0.0,
        "phase8c_mean_cm_dist": -1.0,
        "phase8c_min_cm_dist": -1.0,
        "phase8c_max_cm_dist": -1.0,
        "phase8c_mean_bp_dist": -1.0,
        "phase8c_min_bp_dist": -1.0,
        "phase8c_max_bp_dist": -1.0,
    }

    if pos not in states:
        return dict(empty)

    left = [a for a in anchors if a[0] < pos]
    right = [a for a in anchors if a[0] > pos]

    all_sc = score_context(pos, anchors, panel, map_bp, map_cm, decay_cm, top_k)
    left_sc = score_context(pos, left, panel, map_bp, map_cm, decay_cm, top_k)
    right_sc = score_context(pos, right, panel, map_bp, map_cm, decay_cm, top_k)

    left_pred = 1 if left_sc["plus_prob"] >= left_sc["minus_prob"] and left_sc["confidence"] > 0 else -1
    right_pred = 1 if right_sc["plus_prob"] >= right_sc["minus_prob"] and right_sc["confidence"] > 0 else -1

    two_sided = 1 if left and right else 0
    agree = 1 if two_sided and left_pred == right_pred else 0

    out = dict(empty)
    out["phase8c_candidate_in_panel"] = 1
    out["phase8c_anchor_count"] = len(anchors)
    out["phase8c_left_anchor_count"] = len(left)
    out["phase8c_right_anchor_count"] = len(right)
    out["phase8c_two_sided"] = two_sided

    out["phase8c_copy_plus_prob"] = all_sc["plus_prob"]
    out["phase8c_copy_minus_prob"] = all_sc["minus_prob"]
    out["phase8c_copy_confidence"] = all_sc["confidence"]
    out["phase8c_copy_margin"] = all_sc["margin"]
    out["phase8c_copy_entropy"] = all_sc["entropy"]
    out["phase8c_copy_effective_haps"] = all_sc["effective_haps"]
    out["phase8c_top_plus_score"] = all_sc["top_plus_score"]
    out["phase8c_top_minus_score"] = all_sc["top_minus_score"]
    out["phase8c_top_score_margin"] = all_sc["top_score_margin"]

    out["phase8c_left_plus_prob"] = left_sc["plus_prob"]
    out["phase8c_right_plus_prob"] = right_sc["plus_prob"]
    out["phase8c_left_right_agree"] = agree
    out["phase8c_left_confidence"] = left_sc["confidence"]
    out["phase8c_right_confidence"] = right_sc["confidence"]

    out["phase8c_mean_cm_dist"] = all_sc["mean_cm_dist"]
    out["phase8c_min_cm_dist"] = all_sc["min_cm_dist"]
    out["phase8c_max_cm_dist"] = all_sc["max_cm_dist"]
    out["phase8c_mean_bp_dist"] = all_sc["mean_bp_dist"]
    out["phase8c_min_bp_dist"] = all_sc["min_bp_dist"]
    out["phase8c_max_bp_dist"] = all_sc["max_bp_dist"]

    out["log_phase8c_copy_margin"] = safe_log1p(out["phase8c_copy_margin"])
    out["log_phase8c_copy_effective_haps"] = safe_log1p(out["phase8c_copy_effective_haps"])
    out["log_phase8c_anchor_count"] = safe_log1p(out["phase8c_anchor_count"])

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--genetic-map", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--max-cm", type=float, default=0.25)
    ap.add_argument("--decay-cm", type=float, default=0.05)
    ap.add_argument("--max-anchors-each-side", type=int, default=16)
    ap.add_argument("--top-k", type=int, default=64)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    train = read_tsv(args.training_tsv)
    variants = load_variants(args.variant_json)
    anchors_by_block = load_phase6_anchors(args.local_calls_phase6c_tsv)
    map_bp, map_cm = load_genetic_map(args.genetic_map)

    needed_positions = set()

    candidate_anchors = {}

    for r in train:
        pos = ival(r.get("pos"), None)
        block = str(r.get("block_id", r.get("selected_block_id", "unassigned")))

        if pos is None:
            continue

        if pos in variants:
            needed_positions.add(pos)

        anchors = choose_anchors(
            pos,
            block,
            anchors_by_block,
            map_bp,
            map_cm,
            args.max_cm,
            args.max_anchors_each_side,
        )

        candidate_anchors[pos] = anchors

        for apos, *_ in anchors:
            if apos in variants:
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

    for r in train:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        block = str(rr.get("block_id", rr.get("selected_block_id", "unassigned")))

        if pos is None:
            feats = {}
        else:
            anchors = candidate_anchors.get(pos, [])
            feats = make_features_for_candidate(
                pos,
                block,
                anchors,
                panel,
                map_bp,
                map_cm,
                args.decay_cm,
                args.top_k,
            )

        rr.update(feats)
        out_rows.append(rr)

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

    print(json.dumps({
        "training_rows": len(train),
        "output_rows": len(out_rows),
        "needed_positions": len(needed_positions),
        "panel_positions_needed": panel["positions_needed"],
        "panel_positions_found": panel["positions_found"],
        "panel_allele_mismatch": panel["allele_mismatch"],
        "panel_nonbiallelic": panel["nonbiallelic"],
        "out": args.out_tsv,
    }, indent=2))

if __name__ == "__main__":
    main()
