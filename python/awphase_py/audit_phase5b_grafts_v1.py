#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import Counter, defaultdict
import pysam

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

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def state(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def usable_block_id(row):
    bid = str(row.get("block_id", "")).strip()
    if bid == "" or bid.lower() in {"unassigned", "none", "na", "."}:
        return ""
    return bid

def load_variants(path):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj
    out = {}
    for r in obj:
        try:
            out[int(r["pos"])] = (str(r["ref_allele"]), str(r["alt_allele"]))
        except Exception:
            pass
    return out

def load_truth(path, sample=None):
    vcf = pysam.VariantFile(path)
    samples = list(vcf.header.samples)
    if not samples:
        raise SystemExit(f"No samples in truth VCF: {path}")
    if sample and sample in samples:
        use_sample = sample
    elif sample and sample not in samples and len(samples) == 1:
        use_sample = samples[0]
    elif sample and sample not in samples:
        raise SystemExit(f"Sample {sample} not found. Available={samples}")
    else:
        use_sample = samples[0]

    truth = {}
    for rec in vcf.fetch():
        if rec.alts is None or len(rec.alts) != 1:
            continue
        sd = rec.samples[use_sample]
        gt = sd.get("GT")
        phased = getattr(sd, "phased", False)
        if not phased or gt is None or len(gt) != 2 or None in gt:
            continue
        if tuple(gt) == (0, 1):
            truth[(int(rec.pos), str(rec.ref), str(rec.alts[0]))] = 1
        elif tuple(gt) == (1, 0):
            truth[(int(rec.pos), str(rec.ref), str(rec.alts[0]))] = -1
    return truth, use_sample

def dist_bin(d):
    if d is None:
        return "missing"
    if d <= 1000:
        return "<=1kb"
    if d <= 5000:
        return "1-5kb"
    if d <= 10000:
        return "5-10kb"
    if d <= 25000:
        return "10-25kb"
    if d <= 50000:
        return "25-50kb"
    if d <= 100000:
        return "50-100kb"
    return ">100kb"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tagged-local-calls-tsv", required=True)
    ap.add_argument("--grafted-local-calls-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--truth-sample")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    tagged = {ival(r["pos"]): r for r in read_tsv(args.tagged_local_calls_tsv)}
    grafted = {ival(r["pos"]): r for r in read_tsv(args.grafted_local_calls_tsv)}
    variants = load_variants(args.variant_json)
    truth, truth_sample = load_truth(args.truth_vcf, args.truth_sample)

    rows = []
    bins = defaultdict(lambda: Counter())

    for pos, gr in sorted(grafted.items()):
        tr = tagged.get(pos)
        if tr is None:
            continue

        tagged_state = state(tr)
        grafted_state = state(gr)

        # Newly grafted/fill sites are where tagged was abstain and final has a state.
        if tagged_state != 0 or grafted_state == 0:
            continue

        refalt = variants.get(pos)
        truth_state = ""
        same = ""
        opposite = ""

        if refalt is not None:
            t = truth.get((pos, refalt[0], refalt[1]))
            if t is not None:
                truth_state = t
                same = int(grafted_state == t)
                opposite = int(grafted_state == -t)

        anchor_pos = gr.get("phase5b_anchor_pos", "")
        anchor_state = gr.get("phase5b_anchor_state", "")
        graft_block = gr.get("block_id", "")
        old_block = gr.get("phase5b_old_block_id", "")

        try:
            anchor_dist = abs(pos - int(float(anchor_pos))) if anchor_pos != "" else None
        except Exception:
            anchor_dist = None

        b = dist_bin(anchor_dist)
        bins[b]["n"] += 1
        if same != "":
            bins[b]["truth_comparable"] += 1
            bins[b]["same"] += int(same)
            bins[b]["opposite"] += int(opposite)

        rows.append({
            "pos": pos,
            "tagged_state": tagged_state,
            "grafted_state": grafted_state,
            "truth_state": truth_state,
            "same_orientation_correct": same,
            "opposite_orientation_correct": opposite,
            "anchor_pos": anchor_pos,
            "anchor_state": anchor_state,
            "anchor_dist_bp": "" if anchor_dist is None else anchor_dist,
            "anchor_dist_bin": b,
            "old_block_id": old_block,
            "grafted_block_id": graft_block,
            "phase5b_mode": gr.get("phase5b_mode", ""),
        })

    write_tsv(args.out_tsv, rows)

    summary_bins = {}
    for b, c in bins.items():
        tc = c["truth_comparable"]
        summary_bins[b] = {
            "n": c["n"],
            "truth_comparable": tc,
            "same_correct": c["same"],
            "opposite_correct": c["opposite"],
            "same_accuracy": c["same"] / tc if tc else None,
            "opposite_accuracy": c["opposite"] / tc if tc else None,
        }

    comparable = [r for r in rows if r["truth_state"] != ""]
    same_total = sum(int(r["same_orientation_correct"]) for r in comparable)
    opp_total = sum(int(r["opposite_orientation_correct"]) for r in comparable)

    summary = {
        "truth_sample": truth_sample,
        "grafted_sites": len(rows),
        "truth_comparable_grafted_sites": len(comparable),
        "same_orientation_correct": same_total,
        "opposite_orientation_correct": opp_total,
        "same_orientation_accuracy": same_total / len(comparable) if comparable else None,
        "opposite_orientation_accuracy": opp_total / len(comparable) if comparable else None,
        "distance_bins": summary_bins,
    }

    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
