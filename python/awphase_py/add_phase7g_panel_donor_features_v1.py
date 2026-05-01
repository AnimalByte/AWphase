#!/usr/bin/env python3
import argparse
import csv
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

def safe_log1p(x):
    x = fval(x, 0.0)
    if x < 0:
        x = 0.0
    return math.log1p(x)

def entropy_binary(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))

def load_phase6_anchors(local_calls_tsv):
    rows = read_tsv(local_calls_tsv)
    anchors_by_block = defaultdict(list)
    all_anchors = []

    for r in rows:
        pos = ival(r.get("pos"), None)
        state = ival(r.get("local_phase_state", r.get("phase_state", 0)), 0)
        block = str(r.get("block_id", "unassigned"))

        if pos is None:
            continue

        if state != 0 and block not in {"", "0", "unassigned"}:
            item = (pos, state, block)
            anchors_by_block[block].append(item)
            all_anchors.append(item)

    for b in anchors_by_block:
        anchors_by_block[b].sort(key=lambda x: x[0])
    all_anchors.sort(key=lambda x: x[0])
    return anchors_by_block, all_anchors

def choose_anchors(pos, block, anchors_by_block, max_dist_bp, max_each_side):
    arr = anchors_by_block.get(block, [])
    if not arr:
        return []

    poss = [x[0] for x in arr]
    i = bisect_left(poss, pos)

    chosen = []

    left = []
    j = i - 1
    while j >= 0 and pos - arr[j][0] <= max_dist_bp and len(left) < max_each_side:
        left.append(arr[j])
        j -= 1

    right = []
    j = i
    while j < len(arr) and arr[j][0] - pos <= max_dist_bp and len(right) < max_each_side:
        right.append(arr[j])
        j += 1

    chosen.extend(reversed(left))
    chosen.extend(right)
    return chosen

def load_panel_states(panel_bcf, chrom, start, end, needed_variants):
    """
    Store phased heterozygous panel states at needed positions.
    state_arrays[pos] = (a0, a1, valid_het)
      a0/a1: uint8 allele on hap0/hap1
      valid_het: bool, phased het biallelic only
    """
    vf = pysam.VariantFile(panel_bcf)

    samples = list(vf.header.samples)
    n = len(samples)
    needed_positions = set(needed_variants)

    state_arrays = {}
    found = 0
    mismatch = 0
    nonbiallelic = 0

    # pysam fetch is 0-based half-open.
    for rec in vf.fetch(chrom, max(0, start - 1), end):
        pos = int(rec.pos)
        if pos not in needed_positions:
            continue

        want_ref, want_alt = needed_variants[pos]
        ref = str(rec.ref).upper()
        alts = rec.alts or []

        if len(alts) != 1:
            nonbiallelic += 1
            continue

        alt = str(alts[0]).upper()
        if ref != want_ref or alt != want_alt:
            mismatch += 1
            continue

        a0 = np.zeros(n, dtype=np.uint8)
        a1 = np.zeros(n, dtype=np.uint8)
        valid = np.zeros(n, dtype=bool)

        for idx, sample in enumerate(samples):
            s = rec.samples[sample]
            gt = s.get("GT")
            phased = bool(s.phased)

            if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None:
                continue
            if not phased:
                continue
            if gt[0] not in (0, 1) or gt[1] not in (0, 1):
                continue

            # Need phased heterozygotes for phase-relation evidence.
            if gt[0] + gt[1] != 1:
                continue

            a0[idx] = gt[0]
            a1[idx] = gt[1]
            valid[idx] = True

        state_arrays[pos] = (a0, a1, valid)
        found += 1

    vf.close()

    return {
        "samples": samples,
        "n_samples": n,
        "states": state_arrays,
        "positions_needed": len(needed_positions),
        "positions_found": found,
        "allele_mismatch_records": mismatch,
        "nonbiallelic_records": nonbiallelic,
    }

def donor_features_for_candidate(pos, block, candidate_anchors, panel_states, recomb_scale_bp):
    states = panel_states["states"]
    n_samples = panel_states["n_samples"]

    empty = {
        "donor_candidate_in_panel": 0,
        "donor_anchor_count": 0,
        "donor_panel_usable_anchors": 0,
        "donor_support_plus": 0.0,
        "donor_support_minus": 0.0,
        "donor_total_support": 0.0,
        "donor_margin": 0.0,
        "donor_confidence": 0.0,
        "donor_entropy": 0.0,
        "donor_top1_score": 0.0,
        "donor_top2_score": 0.0,
        "donor_top_margin": 0.0,
        "donor_effective_n": 0.0,
        "donor_left_support_plus": 0.0,
        "donor_left_support_minus": 0.0,
        "donor_right_support_plus": 0.0,
        "donor_right_support_minus": 0.0,
        "donor_left_confidence": 0.0,
        "donor_right_confidence": 0.0,
        "donor_left_right_agree": 0,
        "donor_two_sided": 0,
        "donor_mean_anchor_dist_bp": -1.0,
        "donor_min_anchor_dist_bp": -1.0,
        "donor_max_anchor_dist_bp": -1.0,
        "hmm_like_log_margin": 0.0,
    }

    if pos not in states:
        return empty

    c0, c1, cvalid = states[pos]
    out = dict(empty)
    out["donor_candidate_in_panel"] = 1
    out["donor_anchor_count"] = len(candidate_anchors)

    plus = 0.0
    minus = 0.0
    left_plus = left_minus = right_plus = right_minus = 0.0

    sample_plus = np.zeros(n_samples, dtype=np.float64)
    sample_minus = np.zeros(n_samples, dtype=np.float64)

    used_anchors = 0
    dists = []

    for apos, astate, ablock in candidate_anchors:
        if apos not in states:
            continue

        a0, a1, avalid = states[apos]
        both = cvalid & avalid
        if not np.any(both):
            continue

        same = both & (c0 == a0) & (c1 == a1)
        opposite = both & (c0 == a1) & (c1 == a0)

        if not np.any(same | opposite):
            continue

        dist = abs(pos - apos)
        dists.append(dist)

        # HMM-ish distance decay. This is not a genetic-map HMM yet; it is a BP proxy.
        w = math.exp(-dist / max(recomb_scale_bp, 1.0))

        # Relation candidate-vs-anchor:
        # same panel phase => candidate state agrees with anchor state
        # opposite panel phase => candidate state is opposite anchor state
        if astate == 1:
            plus_mask = same
            minus_mask = opposite
        else:
            plus_mask = opposite
            minus_mask = same

        pcount = float(np.count_nonzero(plus_mask)) * w
        mcount = float(np.count_nonzero(minus_mask)) * w

        plus += pcount
        minus += mcount

        sample_plus[plus_mask] += w
        sample_minus[minus_mask] += w

        if apos < pos:
            left_plus += pcount
            left_minus += mcount
        elif apos > pos:
            right_plus += pcount
            right_minus += mcount

        used_anchors += 1

    total = plus + minus
    if total <= 0.0:
        out["donor_panel_usable_anchors"] = used_anchors
        return out

    major = max(plus, minus)
    minor = min(plus, minus)
    p_major = major / total

    sample_best = np.maximum(sample_plus, sample_minus)
    sample_best = sample_best[sample_best > 0]

    if sample_best.size:
        sorted_scores = np.sort(sample_best)[::-1]
        top1 = float(sorted_scores[0])
        top2 = float(sorted_scores[1]) if sorted_scores.size > 1 else 0.0
        eff_n = float((sample_best.sum() ** 2) / max(np.sum(sample_best ** 2), 1e-12))
    else:
        top1 = top2 = eff_n = 0.0

    left_total = left_plus + left_minus
    right_total = right_plus + right_minus

    left_pred = 0
    right_pred = 0

    if left_total > 0:
        left_pred = 1 if left_plus >= left_minus else -1
    if right_total > 0:
        right_pred = 1 if right_plus >= right_minus else -1

    out["donor_panel_usable_anchors"] = used_anchors
    out["donor_support_plus"] = plus
    out["donor_support_minus"] = minus
    out["donor_total_support"] = total
    out["donor_margin"] = major - minor
    out["donor_confidence"] = p_major
    out["donor_entropy"] = entropy_binary(p_major)
    out["donor_top1_score"] = top1
    out["donor_top2_score"] = top2
    out["donor_top_margin"] = top1 - top2
    out["donor_effective_n"] = eff_n

    out["donor_left_support_plus"] = left_plus
    out["donor_left_support_minus"] = left_minus
    out["donor_right_support_plus"] = right_plus
    out["donor_right_support_minus"] = right_minus
    out["donor_left_confidence"] = max(left_plus, left_minus) / left_total if left_total > 0 else 0.0
    out["donor_right_confidence"] = max(right_plus, right_minus) / right_total if right_total > 0 else 0.0
    out["donor_two_sided"] = 1 if left_total > 0 and right_total > 0 else 0
    out["donor_left_right_agree"] = 1 if left_pred != 0 and right_pred != 0 and left_pred == right_pred else 0

    if dists:
        out["donor_mean_anchor_dist_bp"] = float(sum(dists) / len(dists))
        out["donor_min_anchor_dist_bp"] = float(min(dists))
        out["donor_max_anchor_dist_bp"] = float(max(dists))

    # HMM-like log-odds margin using distance-decayed donor evidence.
    out["hmm_like_log_margin"] = math.log((plus + 1.0) / (minus + 1.0))

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--max-anchor-dist-bp", type=int, default=20000)
    ap.add_argument("--max-anchors-each-side", type=int, default=4)
    ap.add_argument("--recomb-scale-bp", type=float, default=50000.0)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    train = read_tsv(args.training_tsv)
    variants = load_variants(args.variant_json)
    anchors_by_block, all_anchors = load_phase6_anchors(args.local_calls_phase6c_tsv)

    needed_positions = set()

    for r in train:
        pos = ival(r.get("pos"), None)
        if pos is not None and pos in variants:
            needed_positions.add(pos)

        block = str(r.get("block_id", "unassigned"))
        if pos is not None:
            for apos, astate, ablock in choose_anchors(
                pos,
                block,
                anchors_by_block,
                args.max_anchor_dist_bp,
                args.max_anchors_each_side,
            ):
                if apos in variants:
                    needed_positions.add(apos)

    needed_variants = {p: variants[p] for p in needed_positions if p in variants}

    panel = load_panel_states(
        args.panel_bcf,
        args.chrom,
        args.start,
        args.end,
        needed_variants,
    )

    out = []
    for r in train:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        block = str(rr.get("block_id", "unassigned"))

        if pos is None:
            feats = donor_features_for_candidate(
                -1, block, [], panel, args.recomb_scale_bp
            )
        else:
            candidate_anchors = choose_anchors(
                pos,
                block,
                anchors_by_block,
                args.max_anchor_dist_bp,
                args.max_anchors_each_side,
            )
            feats = donor_features_for_candidate(
                pos, block, candidate_anchors, panel, args.recomb_scale_bp
            )

        # Log features
        feats["log_donor_margin"] = safe_log1p(feats["donor_margin"])
        feats["log_donor_total_support"] = safe_log1p(feats["donor_total_support"])
        feats["log_donor_top1_score"] = safe_log1p(feats["donor_top1_score"])
        feats["log_donor_effective_n"] = safe_log1p(feats["donor_effective_n"])
        feats["abs_hmm_like_log_margin"] = abs(feats["hmm_like_log_margin"])

        rr.update(feats)
        out.append(rr)

    fields = []
    seen = set()
    for r in out:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print({
        "training_rows": len(train),
        "output_rows": len(out),
        "needed_positions": len(needed_positions),
        "panel_positions_found": panel["positions_found"],
        "panel_positions_needed": panel["positions_needed"],
        "panel_allele_mismatch_records": panel["allele_mismatch_records"],
        "panel_nonbiallelic_records": panel["nonbiallelic_records"],
        "out": args.out_tsv,
    })

if __name__ == "__main__":
    main()
