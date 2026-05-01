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

def entropy_probs(ps):
    total = float(sum(ps))
    if total <= 0:
        return 0.0
    out = 0.0
    for p in ps:
        q = float(p) / total
        if q > 0:
            out -= q * math.log(q)
    return out

def effective_n(weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w[w > 0]
    if w.size == 0:
        return 0.0
    return float((w.sum() ** 2) / max(float((w * w).sum()), 1e-12))

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
    return (
        np.array([bp[i] for i in order], dtype=np.float64),
        np.array([cm[i] for i in order], dtype=np.float64),
    )

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

def infer_metadata_columns(meta_rows, panel_samples):
    if not meta_rows:
        raise SystemExit("metadata table is empty")

    cols = list(meta_rows[0].keys())
    sample_set = set(panel_samples)

    best_sample_col = None
    best_overlap = -1

    for col in cols:
        vals = set(str(r.get(col, "")).strip() for r in meta_rows)
        overlap = len(vals & sample_set)
        if overlap > best_overlap:
            best_overlap = overlap
            best_sample_col = col

    if best_overlap < 100:
        raise SystemExit(
            f"Could not confidently match metadata samples to panel samples. "
            f"Best sample column={best_sample_col}, overlap={best_overlap}"
        )

    preferred_pop_cols = [
        "genetic_ancestry_group",
        "global_ancestry",
        "gnomad_population",
        "population",
        "pop",
        "superpopulation",
        "super_pop",
        "continental_population",
        "region",
        "group",
        "subpop",
    ]

    lower_to_col = {c.lower(): c for c in cols}
    chosen_pop = None

    for name in preferred_pop_cols:
        if name in lower_to_col and lower_to_col[name] != best_sample_col:
            col = lower_to_col[name]
            vals = [str(r.get(col, "")).strip() for r in meta_rows if str(r.get(col, "")).strip()]
            uniq = set(vals)
            if 2 <= len(uniq) <= 100:
                chosen_pop = col
                break

    if chosen_pop is None:
        candidates = []
        for col in cols:
            if col == best_sample_col:
                continue
            vals = [str(r.get(col, "")).strip() for r in meta_rows if str(r.get(col, "")).strip()]
            uniq = set(vals)
            if 2 <= len(uniq) <= 100:
                name_score = 1 if any(x in col.lower() for x in ["pop", "ancestry", "region", "group"]) else 0
                candidates.append((name_score, -len(uniq), col))
        if not candidates:
            raise SystemExit("Could not infer a usable population/ancestry column from metadata")
        candidates.sort(reverse=True)
        chosen_pop = candidates[0][2]

    return best_sample_col, chosen_pop, best_overlap

def load_sample_pop(meta_tsv, panel_samples):
    rows = read_tsv(meta_tsv)
    sample_col, pop_col, overlap = infer_metadata_columns(rows, panel_samples)

    sample_to_pop = {}
    for r in rows:
        sid = str(r.get(sample_col, "")).strip()
        pop = str(r.get(pop_col, "")).strip()
        if sid and pop:
            sample_to_pop[sid] = pop

    mapped = sum(1 for s in panel_samples if s in sample_to_pop)

    return sample_to_pop, {
        "metadata_sample_col": sample_col,
        "metadata_pop_col": pop_col,
        "metadata_panel_sample_overlap": overlap,
        "panel_samples_mapped_to_pop": mapped,
    }

def ancestry_features(pos, anchors, panel, hap_pop_id, n_pops, pop_hap_counts, map_bp, map_cm, decay_cm):
    states = panel["states"]
    n_haps = panel["n_haps"]

    empty = {
        "phase8e_candidate_in_panel": 0,
        "phase8e_anchor_count": len(anchors),
        "phase8e_usable_anchors": 0,
        "phase8e_pop_n_nonzero": 0,
        "phase8e_pop_top1_prior": 0.0,
        "phase8e_pop_top2_prior": 0.0,
        "phase8e_pop_prior_margin": 0.0,
        "phase8e_pop_prior_entropy": 0.0,
        "phase8e_pop_effective_n": 0.0,
        "phase8e_weighted_plus_support": 0.0,
        "phase8e_weighted_minus_support": 0.0,
        "phase8e_weighted_total_support": 0.0,
        "phase8e_weighted_margin": 0.0,
        "phase8e_weighted_confidence": 0.0,
        "phase8e_weighted_entropy": 0.0,
        "phase8e_weighted_log_margin": 0.0,
        "phase8e_unweighted_plus_support": 0.0,
        "phase8e_unweighted_minus_support": 0.0,
        "phase8e_unweighted_confidence": 0.0,
        "phase8e_weighted_unweighted_delta": 0.0,
    }

    if pos not in states or not anchors:
        return dict(empty)

    cand_alleles, cand_valid = states[pos]

    plus_score = np.zeros(n_haps, dtype=np.float64)
    minus_score = np.zeros(n_haps, dtype=np.float64)
    pop_scores = np.zeros(n_pops, dtype=np.float64)

    qcm = interp_cm(pos, map_bp, map_cm)
    usable = 0

    for apos, astate, *_ in anchors:
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
        any_match = plus_match | minus_match

        plus_score[plus_match] += w
        minus_score[minus_match] += w

        if np.any(any_match):
            pop_scores += np.bincount(
                hap_pop_id[any_match],
                weights=np.full(int(any_match.sum()), w, dtype=np.float64),
                minlength=n_pops,
            )

        usable += 1

    if usable == 0:
        out = dict(empty)
        out["phase8e_candidate_in_panel"] = 1
        return out

    # Normalize ancestry/pop score by available haplotypes per population.
    norm_pop_scores = pop_scores / np.maximum(pop_hap_counts, 1.0)
    total_pop = float(norm_pop_scores.sum())

    if total_pop > 0:
        pop_prior = norm_pop_scores / total_pop
    else:
        pop_prior = np.ones(n_pops, dtype=np.float64) / max(n_pops, 1)

    sorted_prior = np.sort(pop_prior)[::-1]
    top1 = float(sorted_prior[0]) if sorted_prior.size else 0.0
    top2 = float(sorted_prior[1]) if sorted_prior.size > 1 else 0.0

    hap_prior = pop_prior[hap_pop_id]

    plus_mask = cand_valid & (plus_score > 0)
    minus_mask = cand_valid & (minus_score > 0)

    # Candidate state +1:
    # plus-copying hap carries ALT and minus-copying hap carries REF.
    w_plus_support = float(np.sum(plus_score[plus_mask & (cand_alleles == 1)] * hap_prior[plus_mask & (cand_alleles == 1)]))
    w_plus_support += float(np.sum(minus_score[minus_mask & (cand_alleles == 0)] * hap_prior[minus_mask & (cand_alleles == 0)]))

    w_minus_support = float(np.sum(plus_score[plus_mask & (cand_alleles == 0)] * hap_prior[plus_mask & (cand_alleles == 0)]))
    w_minus_support += float(np.sum(minus_score[minus_mask & (cand_alleles == 1)] * hap_prior[minus_mask & (cand_alleles == 1)]))

    u_plus_support = float(np.sum(plus_score[plus_mask & (cand_alleles == 1)]))
    u_plus_support += float(np.sum(minus_score[minus_mask & (cand_alleles == 0)]))

    u_minus_support = float(np.sum(plus_score[plus_mask & (cand_alleles == 0)]))
    u_minus_support += float(np.sum(minus_score[minus_mask & (cand_alleles == 1)]))

    w_total = w_plus_support + w_minus_support
    u_total = u_plus_support + u_minus_support

    if w_total > 0:
        w_conf = max(w_plus_support, w_minus_support) / w_total
        w_margin = abs(w_plus_support - w_minus_support)
        w_log_margin = math.log((w_plus_support + 1.0) / (w_minus_support + 1.0))
        w_entropy = entropy_probs([w_plus_support, w_minus_support])
    else:
        w_conf = w_margin = w_log_margin = w_entropy = 0.0

    if u_total > 0:
        u_conf = max(u_plus_support, u_minus_support) / u_total
    else:
        u_conf = 0.0

    out = dict(empty)
    out["phase8e_candidate_in_panel"] = 1
    out["phase8e_anchor_count"] = len(anchors)
    out["phase8e_usable_anchors"] = usable
    out["phase8e_pop_n_nonzero"] = int(np.sum(pop_prior > 0.001))
    out["phase8e_pop_top1_prior"] = top1
    out["phase8e_pop_top2_prior"] = top2
    out["phase8e_pop_prior_margin"] = top1 - top2
    out["phase8e_pop_prior_entropy"] = entropy_probs(pop_prior)
    out["phase8e_pop_effective_n"] = effective_n(pop_prior)
    out["phase8e_weighted_plus_support"] = w_plus_support
    out["phase8e_weighted_minus_support"] = w_minus_support
    out["phase8e_weighted_total_support"] = w_total
    out["phase8e_weighted_margin"] = w_margin
    out["phase8e_weighted_confidence"] = w_conf
    out["phase8e_weighted_entropy"] = w_entropy
    out["phase8e_weighted_log_margin"] = w_log_margin
    out["phase8e_unweighted_plus_support"] = u_plus_support
    out["phase8e_unweighted_minus_support"] = u_minus_support
    out["phase8e_unweighted_confidence"] = u_conf
    out["phase8e_weighted_unweighted_delta"] = w_conf - u_conf

    out["log_phase8e_weighted_total_support"] = safe_log1p(w_total)
    out["log_phase8e_weighted_margin"] = safe_log1p(w_margin)
    out["log_phase8e_pop_effective_n"] = safe_log1p(out["phase8e_pop_effective_n"])

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-tsv", required=True)
    ap.add_argument("--local-calls-phase6c-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--panel-bcf", required=True)
    ap.add_argument("--panel-meta-tsv", required=True)
    ap.add_argument("--genetic-map", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--max-cm", type=float, default=0.25)
    ap.add_argument("--decay-cm", type=float, default=0.05)
    ap.add_argument("--max-anchors-each-side", type=int, default=16)
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

        candidate_anchors[(pos, block)] = anchors

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

    sample_to_pop, meta_summary = load_sample_pop(args.panel_meta_tsv, panel["samples"])

    pop_labels = sorted(set(sample_to_pop.get(s, "UNKNOWN") for s in panel["samples"]))
    if "UNKNOWN" not in pop_labels:
        pop_labels.append("UNKNOWN")

    pop_to_id = {p: i for i, p in enumerate(pop_labels)}
    n_pops = len(pop_labels)

    hap_pop_id = np.zeros(panel["n_haps"], dtype=np.int64)
    for si, sample in enumerate(panel["samples"]):
        pop = sample_to_pop.get(sample, "UNKNOWN")
        pid = pop_to_id.get(pop, pop_to_id["UNKNOWN"])
        hap_pop_id[2 * si] = pid
        hap_pop_id[2 * si + 1] = pid

    pop_hap_counts = np.bincount(hap_pop_id, minlength=n_pops).astype(np.float64)

    out_rows = []
    for r in train:
        rr = dict(r)
        pos = ival(rr.get("pos"), None)
        block = str(rr.get("block_id", rr.get("selected_block_id", "unassigned")))

        if pos is None:
            feats = {}
        else:
            anchors = candidate_anchors.get((pos, block), [])
            feats = ancestry_features(
                pos,
                anchors,
                panel,
                hap_pop_id,
                n_pops,
                pop_hap_counts,
                map_bp,
                map_cm,
                args.decay_cm,
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

    summary = {
        "training_rows": len(train),
        "output_rows": len(out_rows),
        "needed_positions": len(needed_positions),
        "panel_positions_needed": panel["positions_needed"],
        "panel_positions_found": panel["positions_found"],
        "panel_allele_mismatch": panel["allele_mismatch"],
        "panel_nonbiallelic": panel["nonbiallelic"],
        "n_pop_labels": n_pops,
        "pop_labels": pop_labels,
        "out": args.out_tsv,
    }
    summary.update(meta_summary)

    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
