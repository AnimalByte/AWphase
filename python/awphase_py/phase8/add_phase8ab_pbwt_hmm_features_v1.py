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
    """
    Robust loader for PLINK-style maps:
      chrom  id  cM  bp

    Also tolerates gzipped maps.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    xs = []
    ys = []

    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                cm = float(parts[-2])
                bp = int(float(parts[-1]))
            except Exception:
                continue
            xs.append(bp)
            ys.append(cm)

    if not xs:
        raise SystemExit(f"No map points loaded from {path}")

    order = np.argsort(xs)
    x = np.array([xs[i] for i in order], dtype=np.float64)
    y = np.array([ys[i] for i in order], dtype=np.float64)
    return x, y

def interp_cm(pos, map_bp, map_cm):
    return float(np.interp(float(pos), map_bp, map_cm))

def entropy_binary(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -(p * math.log(p) + (1-p) * math.log(1-p))

def safe_log1p(x):
    x = fval(x, 0.0)
    if x < 0:
        x = 0.0
    return math.log1p(x)

def load_phase6_anchors(path):
    rows = read_tsv(path)
    by_block = defaultdict(list)

    for r in rows:
        pos = ival(r.get("pos"), None)
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))

        if pos is None:
            continue
        if state != 0 and block not in {"", "0", "unassigned"}:
            by_block[block].append((pos, state))

    for b in by_block:
        by_block[b].sort()
    return by_block

def choose_anchors(pos, block, anchors_by_block, map_bp, map_cm, max_cm, max_each_side):
    arr = anchors_by_block.get(block, [])
    if not arr:
        return []

    poss = [x[0] for x in arr]
    cms = [interp_cm(p, map_bp, map_cm) for p in poss]
    qcm = interp_cm(pos, map_bp, map_cm)

    i = bisect_left(poss, pos)
    left = []
    right = []

    j = i - 1
    while j >= 0 and abs(qcm - cms[j]) <= max_cm and len(left) < max_each_side:
        left.append(arr[j])
        j -= 1

    j = i
    while j < len(arr) and abs(cms[j] - qcm) <= max_cm and len(right) < max_each_side:
        right.append(arr[j])
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

            alleles[2*si] = gt[0]
            alleles[2*si + 1] = gt[1]
            valid[2*si] = True
            valid[2*si + 1] = True

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
    return float((w.sum() ** 2) / max(float((w*w).sum()), 1e-12))

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

def phase8_features_for_candidate(pos, block, anchors, panel, map_bp, map_cm, decay_cm, top_k):
    empty = {
        "phase8_candidate_in_panel": 0,
        "phase8_anchor_count": len(anchors),
        "phase8_panel_usable_anchors": 0,
        "phase8_cm_span": 0.0,
        "phase8_pbwt_plus_support": 0.0,
        "phase8_pbwt_minus_support": 0.0,
        "phase8_pbwt_total_support": 0.0,
        "phase8_pbwt_margin": 0.0,
        "phase8_pbwt_confidence": 0.0,
        "phase8_pbwt_entropy": 0.0,
        "phase8_pbwt_top_plus_score": 0.0,
        "phase8_pbwt_top_minus_score": 0.0,
        "phase8_pbwt_top_margin": 0.0,
        "phase8_pbwt_effective_donors": 0.0,
        "phase8_hmm_log_margin": 0.0,
        "phase8_hmm_abs_log_margin": 0.0,
        "phase8_left_support": 0.0,
        "phase8_right_support": 0.0,
        "phase8_two_sided": 0,
        "phase8_left_right_agree": 0,
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
    anchor_cms = []

    usable = 0
    left_support = 0.0
    right_support = 0.0
    left_pred = 0
    right_pred = 0

    for apos, astate in anchors:
        if apos not in states:
            continue

        a_alleles, a_valid = states[apos]
        acm = interp_cm(apos, map_bp, map_cm)
        dcm = abs(qcm - acm)
        w = math.exp(-dcm / max(decay_cm, 1e-9))

        # local plus haplotype expectation:
        # if state=+1, plus carries ALT(1), minus carries REF(0)
        # if state=-1, plus carries REF(0), minus carries ALT(1)
        plus_expected = 1 if astate == 1 else 0
        minus_expected = 1 - plus_expected

        ok = a_valid
        plus_match = ok & (a_alleles == plus_expected)
        minus_match = ok & (a_alleles == minus_expected)

        plus_score[plus_match] += w
        minus_score[minus_match] += w
        plus_valid |= plus_match
        minus_valid |= minus_match

        usable += 1
        anchor_cms.append(acm)

    if usable == 0:
        out = dict(empty)
        out["phase8_candidate_in_panel"] = 1
        return out

    plus_idx, plus_vals = topk_scores(plus_score, plus_valid & cand_valid, top_k)
    minus_idx, minus_vals = topk_scores(minus_score, minus_valid & cand_valid, top_k)

    support_plus = 0.0
    support_minus = 0.0

    # If top plus donors carry ALT at candidate, that supports candidate state +.
    # If top minus donors carry REF at candidate, that also supports candidate state +.
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
    else:
        major = minor = conf = log_margin = 0.0

    all_top_vals = list(plus_vals) + list(minus_vals)
    top_plus = float(plus_vals[0]) if len(plus_vals) else 0.0
    top_minus = float(minus_vals[0]) if len(minus_vals) else 0.0

    out = dict(empty)
    out["phase8_candidate_in_panel"] = 1
    out["phase8_anchor_count"] = len(anchors)
    out["phase8_panel_usable_anchors"] = usable
    out["phase8_cm_span"] = max(anchor_cms) - min(anchor_cms) if anchor_cms else 0.0
    out["phase8_pbwt_plus_support"] = support_plus
    out["phase8_pbwt_minus_support"] = support_minus
    out["phase8_pbwt_total_support"] = total
    out["phase8_pbwt_margin"] = major - minor
    out["phase8_pbwt_confidence"] = conf
    out["phase8_pbwt_entropy"] = entropy_binary(conf)
    out["phase8_pbwt_top_plus_score"] = top_plus
    out["phase8_pbwt_top_minus_score"] = top_minus
    out["phase8_pbwt_top_margin"] = abs(top_plus - top_minus)
    out["phase8_pbwt_effective_donors"] = effective_n(all_top_vals)
    out["phase8_hmm_log_margin"] = log_margin
    out["phase8_hmm_abs_log_margin"] = abs(log_margin)

    out["log_phase8_pbwt_margin"] = safe_log1p(out["phase8_pbwt_margin"])
    out["log_phase8_pbwt_total_support"] = safe_log1p(total)
    out["log_phase8_pbwt_effective_donors"] = safe_log1p(out["phase8_pbwt_effective_donors"])

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
    ap.add_argument("--max-anchors-each-side", type=int, default=12)
    ap.add_argument("--top-k", type=int, default=64)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    train = read_tsv(args.training_tsv)
    variants = load_variants(args.variant_json)
    anchors_by_block = load_phase6_anchors(args.local_calls_phase6c_tsv)
    map_bp, map_cm = load_genetic_map(args.genetic_map)

    needed_positions = set()

    for r in train:
        pos = ival(r.get("pos"), None)
        block = str(r.get("block_id", "unassigned"))
        if pos is None:
            continue
        if pos in variants:
            needed_positions.add(pos)
        anchors = choose_anchors(
            pos, block, anchors_by_block, map_bp, map_cm,
            args.max_cm, args.max_anchors_each_side
        )
        for apos, _ in anchors:
            if apos in variants:
                needed_positions.add(apos)

    needed_variants = {p: variants[p] for p in needed_positions if p in variants}

    panel = load_panel_haps(args.panel_bcf, args.chrom, args.start, args.end, needed_variants)

    out_rows = []
    for r in train:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        block = str(rr.get("block_id", "unassigned"))

        if pos is None:
            feats = {}
        else:
            anchors = choose_anchors(
                pos, block, anchors_by_block, map_bp, map_cm,
                args.max_cm, args.max_anchors_each_side
            )
            feats = phase8_features_for_candidate(
                pos, block, anchors, panel, map_bp, map_cm,
                args.decay_cm, args.top_k
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
