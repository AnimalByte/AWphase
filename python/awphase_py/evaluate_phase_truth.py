#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
from collections import defaultdict

try:
    import pysam
except ImportError as e:
    raise SystemExit(
        "pysam is required. In your hts_extract env, run:\n"
        "  conda install -y -n hts_extract -c conda-forge -c bioconda pysam"
    ) from e


def load_bed_intervals(path):
    intervals = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            intervals[chrom].append((int(start), int(end)))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def interval_span_bp(intervals):
    return sum(
        end - start
        for chrom_intervals in intervals.values()
        for start, end in chrom_intervals
    )


def in_bed(chrom, pos1, intervals):
    zero_pos = pos1 - 1
    for start, end in intervals.get(chrom, []):
        if start <= zero_pos < end:
            return True
        if start > zero_pos:
            return False
    return False


def load_variant_json(path):
    by_pos = {}
    with open(path) as f:
        data = json.load(f)
    for row in data:
        pos = int(row["pos"])
        by_pos[pos] = (str(row["ref_allele"]), str(row["alt_allele"]))
    return by_pos


def truth_state_from_sample(sample_data):
    gt = sample_data.get("GT")
    phased = getattr(sample_data, "phased", False)

    if not phased or gt is None or len(gt) != 2 or None in gt:
        return None, None

    if tuple(gt) == (0, 1):
        return 1, "0|1"
    if tuple(gt) == (1, 0):
        return -1, "1|0"

    return None, None


def load_truth_exact(truth_vcf, truth_bed, sample_name=None):
    intervals = load_bed_intervals(truth_bed)
    vcf = pysam.VariantFile(truth_vcf)

    samples = list(vcf.header.samples)
    if not samples:
        raise SystemExit("No samples found in truth VCF.")
    sample = sample_name or samples[0]
    if sample not in samples:
        raise SystemExit(f"Sample {sample!r} not found in truth VCF. Available: {samples}")

    truth_by_exact = {}
    truth_by_pos = defaultdict(list)

    for chrom, chrom_intervals in intervals.items():
        if chrom not in vcf.header.contigs:
            continue
        fetch_start = min(s for s, _ in chrom_intervals)
        fetch_end = max(e for _, e in chrom_intervals)

        for rec in vcf.fetch(chrom, fetch_start, fetch_end):
            if not in_bed(chrom, rec.pos, intervals):
                continue
            if rec.alts is None or len(rec.alts) != 1:
                continue

            state, gt_str = truth_state_from_sample(rec.samples[sample])
            if state is None:
                continue

            key = (rec.pos, rec.ref, rec.alts[0])
            truth_by_exact[key] = {
                "chrom": chrom,
                "pos": rec.pos,
                "ref": rec.ref,
                "alt": rec.alts[0],
                "truth_state": state,
                "truth_gt": gt_str,
            }
            truth_by_pos[rec.pos].append(key)

    return truth_by_exact, dict(truth_by_pos), sample


def load_pred_rows(pred_tsv, variant_by_pos, phase_column):
    rows = []
    with open(pred_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"pos", "block_id", "confidence", phase_column}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"Missing required TSV columns: {sorted(missing)}; found {reader.fieldnames}"
            )

        for r in reader:
            pos = int(r["pos"])
            ref, alt = variant_by_pos.get(pos, (None, None))
            rows.append(
                {
                    "pos": pos,
                    "block_id": r["block_id"],
                    "pred_state": int(r[phase_column]),
                    "confidence": float(r["confidence"]),
                    "ref": ref,
                    "alt": alt,
                }
            )
    return rows


def compute_n50(lengths):
    return compute_x50(lengths, None)[0]


def compute_x50(lengths, denominator=None):
    lengths = sorted((int(x) for x in lengths if int(x) > 0), reverse=True)
    if not lengths:
        return 0, 0
    total = sum(lengths) if denominator is None else int(denominator)
    if total <= 0:
        return 0, 0
    target = total / 2.0
    running = 0
    for idx, length in enumerate(lengths, start=1):
        running += length
        if running >= target:
            return length, idx
    return 0, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-tsv", required=True)
    ap.add_argument("--variant-json", required=True)
    ap.add_argument("--truth-vcf", required=True)
    ap.add_argument("--truth-bed", required=True)
    ap.add_argument("--sample", default=None)
    ap.add_argument(
        "--phase-column",
        default="stitched_phase_state",
        choices=["stitched_phase_state", "local_phase_state"],
    )
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    variant_by_pos = load_variant_json(args.variant_json)
    truth_intervals = load_bed_intervals(args.truth_bed)
    benchmark_callable_span_bp = interval_span_bp(truth_intervals)
    truth_by_exact, truth_by_pos, truth_sample = load_truth_exact(
        args.truth_vcf, args.truth_bed, args.sample
    )
    pred_rows = load_pred_rows(args.pred_tsv, variant_by_pos, args.phase_column)

    site_rows = []
    matched_rows = []
    positional_overlap = 0
    exact_overlap = 0
    allele_mismatch_positional_overlap = 0

    for row in pred_rows:
        exact_key = None
        truth = None

        if row["ref"] is not None and row["alt"] is not None:
            exact_key = (row["pos"], row["ref"], row["alt"])
            truth = truth_by_exact.get(exact_key)

        pos_truth_keys = truth_by_pos.get(row["pos"], [])
        if pos_truth_keys:
            positional_overlap += 1
            if truth is None:
                allele_mismatch_positional_overlap += 1

        matched_exact = truth is not None
        if matched_exact:
            exact_overlap += 1

        out = {
            "pos": row["pos"],
            "block_id": row["block_id"],
            "pred_state": row["pred_state"],
            "confidence": row["confidence"],
            "ref": row["ref"] or "",
            "alt": row["alt"] or "",
            "matched_exact_truth": int(matched_exact),
            "truth_state": truth["truth_state"] if truth else "",
            "truth_gt": truth["truth_gt"] if truth else "",
            "truth_ref": truth["ref"] if truth else "",
            "truth_alt": truth["alt"] if truth else "",
        }
        site_rows.append(out)

        if matched_exact:
            matched_rows.append(
                {
                    "pos": row["pos"],
                    "block_id": row["block_id"],
                    "pred_state": row["pred_state"],
                    "confidence": row["confidence"],
                    "truth_state": truth["truth_state"],
                    "truth_gt": truth["truth_gt"],
                }
            )

    phased_matched_rows = [
        r for r in matched_rows
        if r["pred_state"] != 0 and r["block_id"] != "unassigned"
    ]

    # Hamming per predicted block, allowing whole-block global flip
    by_block = defaultdict(list)
    for r in phased_matched_rows:
        by_block[r["block_id"]].append(r)

    hamming_errors = 0
    hamming_denom = 0
    switch_errors = 0
    switch_denom = 0
    block_metric_rows = []
    switch_corrected_block_spans = []

    for block_id, rows in by_block.items():
        rows = sorted(rows, key=lambda x: x["pos"])

        as_is = sum(1 for r in rows if r["pred_state"] != r["truth_state"])
        flipped = sum(1 for r in rows if (-r["pred_state"]) != r["truth_state"])

        if flipped < as_is:
            chosen = "flip"
            block_hamming_errors = flipped
        else:
            chosen = "as_is"
            block_hamming_errors = as_is

        hamming_errors += block_hamming_errors
        hamming_denom += len(rows)

        block_switch_errors = 0
        segment_start = rows[0]["pos"]
        for a, b in zip(rows, rows[1:]):
            pred_rel = a["pred_state"] * b["pred_state"]
            truth_rel = a["truth_state"] * b["truth_state"]
            switch_denom += 1
            if pred_rel != truth_rel:
                switch_errors += 1
                block_switch_errors += 1
                switch_corrected_block_spans.append(a["pos"] - segment_start + 1)
                segment_start = b["pos"]
        switch_corrected_block_spans.append(rows[-1]["pos"] - segment_start + 1)

        block_metric_rows.append(
            {
                "block_id": block_id,
                "n_matched_phased_sites": len(rows),
                "start_pos": rows[0]["pos"],
                "end_pos": rows[-1]["pos"],
                "span_bp": rows[-1]["pos"] - rows[0]["pos"] + 1,
                "hamming_errors": block_hamming_errors,
                "hamming_rate": (
                    block_hamming_errors / len(rows) if rows else 0.0
                ),
                "switch_errors": block_switch_errors,
                "switch_corrected_segments": block_switch_errors + 1,
                "chosen_block_polarity": chosen,
            }
        )

    # Raw predicted block N50 from chosen phase column
    pred_block_positions = defaultdict(list)
    for r in pred_rows:
        if r["pred_state"] != 0 and r["block_id"] != "unassigned":
            pred_block_positions[r["block_id"]].append(r["pos"])

    pred_block_spans = []
    for block_id, positions in pred_block_positions.items():
        positions.sort()
        pred_block_spans.append(positions[-1] - positions[0] + 1)

    raw_block_n50_bp, raw_block_l50 = compute_x50(pred_block_spans)
    raw_block_ng50_bp, raw_block_lg50 = compute_x50(
        pred_block_spans, benchmark_callable_span_bp
    )
    switch_corrected_block_n50_bp, switch_corrected_block_l50 = compute_x50(
        switch_corrected_block_spans
    )
    switch_corrected_block_ng50_bp, switch_corrected_block_lg50 = compute_x50(
        switch_corrected_block_spans, benchmark_callable_span_bp
    )

    n_pred_nonzero = sum(1 for r in pred_rows if r["pred_state"] != 0)
    n_exact_overlap_phased = sum(1 for r in matched_rows if r["pred_state"] != 0)

    metrics = {
        "truth_sample": truth_sample,
        "phase_column": args.phase_column,
        "benchmark_callable_span_bp": benchmark_callable_span_bp,
        "n_truth_het_sites_in_bed": len(truth_by_exact),
        "n_pred_rows_total": len(pred_rows),
        "n_pred_sites_nonzero": n_pred_nonzero,
        "n_positional_overlaps": positional_overlap,
        "n_exact_overlaps": exact_overlap,
        "n_allele_mismatch_positional_overlaps": allele_mismatch_positional_overlap,
        "n_exact_overlap_sites_phased": n_exact_overlap_phased,
        "pct_truth_sites_phased": (
            n_exact_overlap_phased / len(truth_by_exact) if truth_by_exact else 0.0
        ),
        "pct_exact_overlap_sites_phased": (
            n_exact_overlap_phased / exact_overlap if exact_overlap else 0.0
        ),
        "hamming_errors": hamming_errors,
        "hamming_denominator": hamming_denom,
        "hamming_error_rate": (
            hamming_errors / hamming_denom if hamming_denom else 0.0
        ),
        "phased_site_accuracy_pct": (
            100.0 * (1.0 - (hamming_errors / hamming_denom)) if hamming_denom else 0.0
        ),
        "truth_correct_pct": (
            100.0 * ((hamming_denom - hamming_errors) / len(truth_by_exact))
            if truth_by_exact else 0.0
        ),
        "switch_errors": switch_errors,
        "switch_denominator": switch_denom,
        "switch_error_rate": (
            switch_errors / switch_denom if switch_denom else 0.0
        ),
        "n_pred_blocks_nonzero": len(pred_block_positions),
        "raw_block_n50_bp": raw_block_n50_bp,
        "raw_block_l50": raw_block_l50,
        "raw_block_ng50_bp": raw_block_ng50_bp,
        "raw_block_lg50": raw_block_lg50,
        "max_block_span_bp": max(pred_block_spans) if pred_block_spans else 0,
        "median_block_span_bp": (
            sorted(pred_block_spans)[len(pred_block_spans) // 2]
            if pred_block_spans
            else 0
        ),
        "switch_corrected_block_n50_bp": switch_corrected_block_n50_bp,
        "switch_corrected_block_l50": switch_corrected_block_l50,
        "switch_corrected_block_ng50_bp": switch_corrected_block_ng50_bp,
        "switch_corrected_block_lg50": switch_corrected_block_lg50,
        "max_switch_corrected_block_span_bp": (
            max(switch_corrected_block_spans) if switch_corrected_block_spans else 0
        ),
        "median_switch_corrected_block_span_bp": (
            sorted(switch_corrected_block_spans)[len(switch_corrected_block_spans) // 2]
            if switch_corrected_block_spans
            else 0
        ),
        "n_blocks_with_truth_overlap": len(by_block),
    }

    metrics_path = args.out_prefix + ".metrics.json"
    site_path = args.out_prefix + ".site_comparison.tsv"
    block_path = args.out_prefix + ".block_metrics.tsv"

    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    with open(site_path, "w", newline="") as f:
        fieldnames = [
            "pos",
            "block_id",
            "pred_state",
            "confidence",
            "ref",
            "alt",
            "matched_exact_truth",
            "truth_state",
            "truth_gt",
            "truth_ref",
            "truth_alt",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    with open(block_path, "w", newline="") as f:
        fieldnames = [
            "block_id",
            "n_matched_phased_sites",
            "start_pos",
            "end_pos",
            "span_bp",
            "hamming_errors",
            "hamming_rate",
            "switch_errors",
            "switch_corrected_segments",
            "chosen_block_polarity",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(block_metric_rows, key=lambda x: (x["start_pos"], x["block_id"])))

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
