#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
from collections import defaultdict
import pysam

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def ival(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def state(row):
    return ival(row.get("local_phase_state", row.get("phase_state", 0)), 0)

def load_variant_positions(path, chrom):
    obj = json.load(open(path))
    obj = obj.get("variants", obj) if isinstance(obj, dict) else obj
    out = set()
    for r in obj:
        try:
            pos = int(r["pos"])
        except Exception:
            continue
        out.add(pos)
    return out

def load_call_positions(path):
    all_pos = set()
    nonzero_pos = set()
    for r in read_tsv(path):
        pos = ival(r.get("pos"))
        if pos <= 0:
            continue
        all_pos.add(pos)
        if state(r) != 0:
            nonzero_pos.add(pos)
    return all_pos, nonzero_pos

def load_obs_positions(path):
    if not path or not Path(path).exists():
        return set(), defaultdict(int)

    positions = set()
    counts = defaultdict(int)

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            pos = ival(r.get("site_pos", r.get("pos", 0)))
            if pos <= 0:
                continue
            positions.add(pos)
            counts[pos] += 1
    return positions, counts

def load_bed_intervals(path, chrom):
    intervals = []
    if not path or not Path(path).exists():
        return intervals

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            c = parts[0]
            if c != chrom:
                continue
            # BED is 0-based start, half-open. Convert to 1-based-ish inclusive for simple filtering.
            start = int(parts[1]) + 1
            end = int(parts[2])
            intervals.append((start, end))
    return intervals

def count_bed_overlap_bp(start, end, intervals):
    total = 0
    for a, b in intervals:
        if b < start:
            continue
        if a > end:
            break
        total += max(0, min(end, b) - max(start, a) + 1)
    return total

def load_truth_phased_positions(vcf_path, chrom, sample=None):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    if not samples:
        return set(), None

    if sample and sample in samples:
        use_sample = sample
    else:
        use_sample = samples[0]

    positions = set()

    # Try region fetch first, then full scan fallback.
    try:
        records = vcf.fetch(chrom)
    except Exception:
        records = vcf.fetch()

    for rec in records:
        if rec.contig != chrom:
            continue
        if rec.alts is None or len(rec.alts) != 1:
            continue

        sd = rec.samples[use_sample]
        gt = sd.get("GT")
        phased = getattr(sd, "phased", False)

        if not phased or gt is None or len(gt) != 2 or None in gt:
            continue

        # Only het phased.
        if tuple(gt) in {(0, 1), (1, 0)}:
            positions.add(int(rec.pos))

    return positions, use_sample

def bin_positions(positions, start, end, window_bp, step_bp):
    positions = sorted(p for p in positions if start <= p <= end)
    out = []
    i = start
    pos_idx = 0

    while i + window_bp - 1 <= end:
        w_start = i
        w_end = i + window_bp - 1

        # Count with simple scan.
        while pos_idx < len(positions) and positions[pos_idx] < w_start:
            pos_idx += 1
        j = pos_idx
        n = 0
        while j < len(positions) and positions[j] <= w_end:
            n += 1
            j += 1

        out.append((w_start, w_end, n))
        i += step_bp

    return out

def count_in_window(pos_set, start, end):
    return sum(1 for p in pos_set if start <= p <= end)

def sum_obs_in_window(obs_counts, start, end):
    return sum(v for p, v in obs_counts.items() if start <= p <= end)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--truth-bed")
    ap.add_argument("--tagged-local-calls-tsv", required=True)
    ap.add_argument("--fill-local-calls-tsv")
    ap.add_argument("--obs-tsv")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)

    ap.add_argument("--chrom-start", type=int, default=1)
    ap.add_argument("--chrom-end", type=int, required=True)
    ap.add_argument("--window-bp", type=int, default=5000000)
    ap.add_argument("--step-bp", type=int, default=1000000)
    ap.add_argument("--min-start-bp", type=int, default=10000000)
    ap.add_argument("--truth-sample", default=None)
    args = ap.parse_args()

    chrom = args.chrom

    variant_pos = load_variant_positions(args.variant_json, chrom)
    truth_pos, truth_sample = load_truth_phased_positions(args.truth_vcf, chrom, args.truth_sample)
    tagged_all, tagged_nonzero = load_call_positions(args.tagged_local_calls_tsv)

    fill_all, fill_nonzero = set(), set()
    if args.fill_local_calls_tsv and Path(args.fill_local_calls_tsv).exists():
        fill_all, fill_nonzero = load_call_positions(args.fill_local_calls_tsv)

    obs_pos, obs_counts = load_obs_positions(args.obs_tsv)
    bed_intervals = load_bed_intervals(args.truth_bed, chrom) if args.truth_bed else []

    starts = range(max(args.chrom_start, args.min_start_bp), args.chrom_end - args.window_bp + 2, args.step_bp)

    rows = []
    for start in starts:
        end = start + args.window_bp - 1

        bed_bp = count_bed_overlap_bp(start, end, bed_intervals) if bed_intervals else 0
        bed_frac = bed_bp / args.window_bp if bed_intervals else None

        n_variant = count_in_window(variant_pos, start, end)
        n_truth = count_in_window(truth_pos, start, end)
        n_tagged_all = count_in_window(tagged_all, start, end)
        n_tagged_nonzero = count_in_window(tagged_nonzero, start, end)
        n_fill_all = count_in_window(fill_all, start, end)
        n_fill_nonzero = count_in_window(fill_nonzero, start, end)
        n_obs_sites = count_in_window(obs_pos, start, end)
        n_obs = sum_obs_in_window(obs_counts, start, end)

        # Heuristic score: high truth + calls + obs + BED coverage.
        score = (
            2.0 * n_truth
            + 1.5 * n_tagged_nonzero
            + 1.0 * n_fill_nonzero
            + 0.02 * n_obs
            + 500.0 * (bed_frac if bed_frac is not None else 0.0)
        )

        rows.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "region": f"{chrom}:{start}-{end}",
            "window_bp": args.window_bp,
            "bed_overlap_bp": bed_bp,
            "bed_overlap_frac": "" if bed_frac is None else f"{bed_frac:.6f}",
            "variant_json_sites": n_variant,
            "truth_phased_het_sites": n_truth,
            "tagged_sites": n_tagged_all,
            "tagged_nonzero_sites": n_tagged_nonzero,
            "fill_sites": n_fill_all,
            "fill_nonzero_sites": n_fill_nonzero,
            "obs_sites": n_obs_sites,
            "obs_records": n_obs,
            "score": f"{score:.3f}",
        })

    rows.sort(key=lambda r: float(r["score"]), reverse=True)

    fields = [
        "chrom", "start", "end", "region", "window_bp",
        "bed_overlap_bp", "bed_overlap_frac",
        "variant_json_sites", "truth_phased_het_sites",
        "tagged_sites", "tagged_nonzero_sites",
        "fill_sites", "fill_nonzero_sites",
        "obs_sites", "obs_records", "score",
    ]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "chrom": chrom,
        "truth_sample": truth_sample,
        "truth_phased_het_sites_total": len(truth_pos),
        "variant_json_sites_total": len(variant_pos),
        "tagged_nonzero_sites_total": len(tagged_nonzero),
        "fill_nonzero_sites_total": len(fill_nonzero),
        "obs_sites_total": len(obs_pos),
        "windows_scored": len(rows),
        "top_region": rows[0]["region"] if rows else None,
        "top_score": rows[0]["score"] if rows else None,
    }

    Path(args.out_json).write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    print()
    print("top windows:")
    for r in rows[:10]:
        print("\t".join(str(r[k]) for k in fields))

if __name__ == "__main__":
    main()
