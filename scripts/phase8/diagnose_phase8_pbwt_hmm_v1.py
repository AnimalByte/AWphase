#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)


METRIC_KEYS = [
    "n_truth_het_sites_in_bed",
    "n_pred_sites_nonzero",
    "n_exact_overlap_sites_phased",
    "hamming_errors",
    "hamming_denominator",
    "hamming_error_rate",
    "switch_errors",
    "switch_denominator",
    "switch_error_rate",
    "truth_correct_pct",
    "phased_site_accuracy_pct",
    "raw_block_n50_bp",
]

SUMMARY_KEYS = [
    "backbone_anchors",
    "candidate_unphased_positions",
    "candidate_block_pairs",
    "panel_positions_needed",
    "panel_positions_found",
    "panel_positions_missing",
    "panel_allele_mismatch_records",
    "panel_nonbiallelic_records",
    "pbwt_index_positions",
    "accepted_sites",
]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path):
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with open(path) as fh:
        return json.load(fh)


def fval(value, default=0.0):
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(value)
    except Exception:
        return default


def ival(value, default=0):
    try:
        if value is None or str(value).strip() == "":
            return default
        return int(float(value))
    except Exception:
        return default


def sval(value):
    return "" if value is None else str(value).strip()


def load_manifest(path, split):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row = {k: sval(v) for k, v in row.items()}
            if split != "all" and row.get("split") != split:
                continue
            rows.append(row)
    return rows


def result_root(label):
    return Path("results/phase7a_windows") / f"{label}_illumina30x_phase7a"


def orientation_by_block(site_rows):
    by_block = defaultdict(list)
    for row in site_rows:
        block = sval(row.get("block_id"))
        pred = ival(row.get("pred_state"), 0)
        truth = ival(row.get("truth_state"), 0)
        matched = sval(row.get("matched_exact_truth")).lower()
        if not block or block in {"unassigned", "none", "na", "."}:
            continue
        if matched in {"0", "false"} or pred == 0 or truth == 0:
            continue
        by_block[block].append((pred, truth))

    orientation = {}
    counts = {}
    for block, pairs in by_block.items():
        as_is = sum(1 for pred, truth in pairs if pred != truth)
        flipped = sum(1 for pred, truth in pairs if -pred != truth)
        orient = -1 if flipped < as_is else 1
        same = sum(1 for pred, truth in pairs if pred == truth)
        opposite = sum(1 for pred, truth in pairs if -pred == truth)
        orientation[block] = orient
        counts[block] = {
            "same": same,
            "opposite": opposite,
            "orientation": orient,
            "orientation_margin": abs(same - opposite),
            "n_orient_sites": len(pairs),
        }
    return orientation, counts


def oriented_correct_by_pos(site_rows):
    orientation, counts = orientation_by_block(site_rows)
    out = {}
    for row in site_rows:
        pos = ival(row.get("pos"), -1)
        block = sval(row.get("block_id"))
        pred = ival(row.get("pred_state"), 0)
        truth = ival(row.get("truth_state"), 0)
        matched = sval(row.get("matched_exact_truth")).lower()
        if matched in {"0", "false"} or pred == 0 or truth == 0:
            continue
        if not block or block in {"unassigned", "none", "na", "."}:
            continue
        orient = orientation.get(block)
        if orient is None:
            continue
        c = counts.get(block, {})
        out[pos] = {
            "truth_state": truth,
            "block_orientation": orient,
            "oriented_pred_state": pred * orient,
            "correct_after_block_orientation": int(pred * orient == truth),
            "orientation_sites": c.get("n_orient_sites", 0),
            "orientation_margin": c.get("orientation_margin", 0),
        }
    return out


def bin_value(feature, value):
    x = fval(value, 0.0)
    if feature == "anchors":
        if x < 2:
            return "lt2"
        if x < 4:
            return "2_3"
        if x < 8:
            return "4_7"
        if x < 16:
            return "8_15"
        return "ge16"
    if feature == "donors":
        if x < 12:
            return "lt12"
        if x < 32:
            return "12_31"
        if x < 64:
            return "32_63"
        return "ge64"
    if feature == "confidence":
        if x < 0.82:
            return "lt0.82"
        if x < 0.90:
            return "0.82_0.90"
        if x < 0.98:
            return "0.90_0.98"
        return "ge0.98"
    if feature == "margin":
        if x < 1.0:
            return "lt1"
        if x < 10.0:
            return "1_10"
        if x < 100.0:
            return "10_100"
        return "ge100"
    return "all"


def summarize_window(manifest_row):
    label = manifest_row["label"]
    split = manifest_row["split"]
    root = result_root(label)
    metrics = load_json(root / "truth_eval_phase8pbwt_hmm.metrics.json")
    summary = load_json(root / "phase8_pbwt_hmm" / "summary.pbwt_hmm.json")

    local_path = root / "local_calls.phase8pbwt_hmm.tsv"
    cand_path = root / "phase8_pbwt_hmm" / "candidates.pbwt_hmm.tsv"
    site_path = root / "truth_eval_phase8pbwt_hmm.site_comparison.tsv"

    if not local_path.exists() or not cand_path.exists() or not site_path.exists():
        return None, [], []

    local_rows = read_tsv(local_path)
    cand_rows = read_tsv(cand_path)
    site_rows = read_tsv(site_path)
    correct_by_pos = oriented_correct_by_pos(site_rows)

    candidate_positions = {ival(r.get("pos"), -1) for r in cand_rows}
    accepted_candidates = [r for r in cand_rows if ival(r.get("accepted"), 0) == 1]

    filled_rows = []
    feature_counts = defaultdict(Counter)
    counts = Counter()

    for row in local_rows:
        if ival(row.get("phase8_pbwt_hmm_filled"), 0) != 1:
            continue
        pos = ival(row.get("pos"), -1)
        method = sval(row.get("phase8_pbwt_hmm_method")) or "unknown"
        detail = {
            "split": split,
            "label": label,
            "pos": pos,
            "block_id": sval(row.get("block_id")),
            "pred_state": ival(row.get("local_phase_state"), 0),
            "method": method,
            "confidence": fval(row.get("phase8_pbwt_hmm_confidence"), 0.0),
            "margin": fval(row.get("phase8_pbwt_hmm_margin"), 0.0),
            "donors": ival(row.get("phase8_pbwt_hmm_donors"), 0),
            "anchors": ival(row.get("phase8_pbwt_hmm_anchors"), 0),
            "comparable": 0,
            "correct_after_block_orientation": "",
            "truth_state": "",
            "block_orientation": "",
            "orientation_sites": "",
            "orientation_margin": "",
        }
        corr = correct_by_pos.get(pos)
        if corr is not None:
            detail.update(corr)
            detail["comparable"] = 1
            if detail["correct_after_block_orientation"]:
                counts[f"filled_correct_{method}"] += 1
            else:
                counts[f"filled_wrong_{method}"] += 1
                if detail["anchors"] < 4:
                    counts["filled_wrong_anchor_lt4"] += 1
                if detail["donors"] < 32:
                    counts["filled_wrong_donor_lt32"] += 1
                if detail["confidence"] < 0.98:
                    counts["filled_wrong_confidence_lt0.98"] += 1
            counts["filled_comparable"] += 1
        else:
            counts["filled_not_comparable"] += 1
        counts["filled_total"] += 1
        counts[f"filled_{method}"] += 1
        filled_rows.append(detail)

        for feature in ["method", "anchors", "donors", "confidence", "margin"]:
            bucket = method if feature == "method" else bin_value(feature, detail[feature])
            key = (split, method, feature, bucket)
            feature_counts[key]["total"] += 1
            if detail["comparable"]:
                feature_counts[key]["comparable"] += 1
                if detail["correct_after_block_orientation"]:
                    feature_counts[key]["correct"] += 1
                else:
                    feature_counts[key]["wrong"] += 1

    row = {
        "split": split,
        "label": label,
        "chrom": manifest_row.get("chrom", ""),
        "start": manifest_row.get("start", ""),
        "end": manifest_row.get("end", ""),
        "candidate_rows": len(cand_rows),
        "candidate_positions": len(candidate_positions),
        "accepted_candidate_rows": len(accepted_candidates),
    }
    for key in METRIC_KEYS:
        row[key] = metrics.get(key, "")
    for key in SUMMARY_KEYS:
        row[key] = summary.get(key, "")
    for key in [
        "filled_total",
        "filled_comparable",
        "filled_not_comparable",
        "filled_pbwt",
        "filled_pbwt_hmm",
        "filled_correct_pbwt",
        "filled_wrong_pbwt",
        "filled_correct_pbwt_hmm",
        "filled_wrong_pbwt_hmm",
        "filled_wrong_anchor_lt4",
        "filled_wrong_donor_lt32",
        "filled_wrong_confidence_lt0.98",
    ]:
        row[key] = counts.get(key, 0)
    comp = counts.get("filled_comparable", 0)
    correct = counts.get("filled_correct_pbwt", 0) + counts.get("filled_correct_pbwt_hmm", 0)
    row["filled_oriented_accuracy"] = correct / comp if comp else ""

    bin_rows = []
    for (split_name, method, feature, bucket), c in feature_counts.items():
        comparable = c.get("comparable", 0)
        bin_rows.append(
            {
                "split": split_name,
                "label": label,
                "method": method,
                "feature": feature,
                "bucket": bucket,
                "total": c.get("total", 0),
                "comparable": comparable,
                "correct": c.get("correct", 0),
                "wrong": c.get("wrong", 0),
                "accuracy": c.get("correct", 0) / comparable if comparable else "",
            }
        )

    return row, filled_rows, bin_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        default="results/phase8f_manifests/phase8f_windows.window_beds.tsv",
    )
    ap.add_argument("--split", default="all", choices=["all", "train", "holdout_reserved"])
    ap.add_argument("--out-dir", default="results/phase8/pbwt_hmm_diagnostics")
    args = ap.parse_args()

    manifest_rows = load_manifest(args.manifest, args.split)
    window_rows = []
    filled_rows = []
    bin_rows = []
    for manifest_row in manifest_rows:
        window, filled, bins = summarize_window(manifest_row)
        if window is None:
            print(f"skip missing Phase8 outputs: {manifest_row.get('label')}", file=sys.stderr)
            continue
        window_rows.append(window)
        filled_rows.extend(filled)
        bin_rows.extend(bins)

    out_dir = Path(args.out_dir)
    window_fields = [
        "split",
        "label",
        "chrom",
        "start",
        "end",
    ] + METRIC_KEYS + SUMMARY_KEYS + [
        "candidate_rows",
        "candidate_positions",
        "accepted_candidate_rows",
        "filled_total",
        "filled_comparable",
        "filled_not_comparable",
        "filled_pbwt",
        "filled_pbwt_hmm",
        "filled_correct_pbwt",
        "filled_wrong_pbwt",
        "filled_correct_pbwt_hmm",
        "filled_wrong_pbwt_hmm",
        "filled_wrong_anchor_lt4",
        "filled_wrong_donor_lt32",
        "filled_wrong_confidence_lt0.98",
        "filled_oriented_accuracy",
    ]
    filled_fields = [
        "split",
        "label",
        "pos",
        "block_id",
        "pred_state",
        "method",
        "confidence",
        "margin",
        "donors",
        "anchors",
        "comparable",
        "correct_after_block_orientation",
        "truth_state",
        "block_orientation",
        "oriented_pred_state",
        "orientation_sites",
        "orientation_margin",
    ]
    bin_fields = [
        "split",
        "label",
        "method",
        "feature",
        "bucket",
        "total",
        "comparable",
        "correct",
        "wrong",
        "accuracy",
    ]

    window_rows.sort(key=lambda r: (fval(r.get("truth_correct_pct"), 0.0), -fval(r.get("hamming_error_rate"), 0.0)))
    write_tsv(out_dir / f"phase8_pbwt_hmm.{args.split}.window_diagnostics.tsv", window_rows, window_fields)
    write_tsv(out_dir / f"phase8_pbwt_hmm.{args.split}.filled_site_diagnostics.tsv", filled_rows, filled_fields)
    write_tsv(out_dir / f"phase8_pbwt_hmm.{args.split}.feature_bins.tsv", bin_rows, bin_fields)
    write_tsv(out_dir / f"phase8_pbwt_hmm.{args.split}.weak_windows.tsv", window_rows[:10], window_fields)

    print(f"wrote {out_dir / f'phase8_pbwt_hmm.{args.split}.window_diagnostics.tsv'}")
    print(f"wrote {out_dir / f'phase8_pbwt_hmm.{args.split}.filled_site_diagnostics.tsv'}")
    print(f"wrote {out_dir / f'phase8_pbwt_hmm.{args.split}.feature_bins.tsv'}")


if __name__ == "__main__":
    main()
