#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from collections import defaultdict

DEFAULT_METHOD_SETS = {
    "common_with_illumina_whatshap": [
        "AWPhase_Phase6C_read_backbone",
        "AWPhase_Phase7A_hand_panel_fill",
        "AWPhase_Phase8F_t0.80",
        "AWPhase_Phase9A_bridge_t0.995_on_Phase8F_t0.80",
        "WhatsHap_Illumina30x",
    ],
    "common_with_hifi_whatshap": [
        "AWPhase_Phase6C_read_backbone",
        "AWPhase_Phase7A_hand_panel_fill",
        "AWPhase_Phase8F_t0.80",
        "AWPhase_Phase9A_bridge_t0.995_on_Phase8F_t0.80",
        "WhatsHap_GIAB_CCS_HiFi30x",
    ],
    "common_all_methods": [
        "AWPhase_Phase6C_read_backbone",
        "AWPhase_Phase7A_hand_panel_fill",
        "AWPhase_Phase8F_t0.80",
        "AWPhase_Phase9A_bridge_t0.995_on_Phase8F_t0.80",
        "WhatsHap_Illumina30x",
        "WhatsHap_GIAB_CCS_HiFi30x",
    ],
}

OUT_FIELDS = [
    "comparison_set",
    "method",
    "n_windows",
    "windows",
    "total_truth_sites",
    "total_pred_sites_nonzero",
    "total_exact_overlap_sites_phased",
    "total_hamming_errors",
    "weighted_hamming_error_rate",
    "total_switch_errors",
    "weighted_switch_error_rate",
    "weighted_phased_site_accuracy_pct",
    "weighted_truth_correct_pct",
    "mean_truth_correct_pct",
    "mean_hamming_error_rate",
    "mean_switch_error_rate",
    "mean_raw_block_n50_bp",
    "mean_median_block_span_bp",
    "mean_max_block_span_bp",
]

REQUIRED_COLUMNS = [
    "label",
    "method",
    "n_truth_het_sites_in_bed",
    "n_pred_sites_nonzero",
    "n_exact_overlap_sites_phased",
    "hamming_errors",
    "hamming_denominator",
    "hamming_error_rate",
    "switch_errors",
    "switch_denominator",
    "switch_error_rate",
    "phased_site_accuracy_pct",
    "truth_correct_pct",
    "raw_block_n50_bp",
    "median_block_span_bp",
    "max_block_span_bp",
]


def f(x):
    try:
        if x is None or str(x).strip() == "":
            return 0.0
        return float(x)
    except Exception:
        return 0.0


def read_tsv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"No rows found in {path}")

    missing = [c for c in REQUIRED_COLUMNS if c not in rows[0]]
    if missing:
        raise SystemExit(f"Input TSV is missing required columns: {missing}")

    return rows


def summarize(rows, methods, comparison_set):
    labels_by_method = defaultdict(set)

    for r in rows:
        if r["method"] in methods:
            labels_by_method[r["method"]].add(r["label"])

    common = None
    for method in methods:
        s = labels_by_method.get(method, set())
        if common is None:
            common = set(s)
        else:
            common &= s

    common = common or set()

    if not common:
        print(f"WARNING: no common windows for {comparison_set}")
        for method in methods:
            print(f"  {method}: {len(labels_by_method.get(method, set()))} windows")
        return [], common

    by_method = defaultdict(list)
    for r in rows:
        if r["method"] in methods and r["label"] in common:
            by_method[r["method"]].append(r)

    out_rows = []

    for method in methods:
        rs = by_method.get(method, [])
        if not rs:
            continue

        total_truth = sum(f(r["n_truth_het_sites_in_bed"]) for r in rs)
        total_pred = sum(f(r["n_pred_sites_nonzero"]) for r in rs)
        total_phased = sum(f(r["n_exact_overlap_sites_phased"]) for r in rs)

        total_hamming = sum(f(r["hamming_errors"]) for r in rs)
        total_hamming_den = sum(f(r["hamming_denominator"]) for r in rs)

        total_switch = sum(f(r["switch_errors"]) for r in rs)
        total_switch_den = sum(f(r["switch_denominator"]) for r in rs)

        weighted_hamming_error_rate = (
            total_hamming / total_hamming_den if total_hamming_den else ""
        )
        weighted_switch_error_rate = (
            total_switch / total_switch_den if total_switch_den else ""
        )
        weighted_phased_site_accuracy_pct = (
            100.0 * (1.0 - total_hamming / total_hamming_den)
            if total_hamming_den
            else ""
        )
        weighted_truth_correct_pct = (
            100.0 * (total_phased - total_hamming) / total_truth
            if total_truth
            else ""
        )

        out_rows.append({
            "comparison_set": comparison_set,
            "method": method,
            "n_windows": len(rs),
            "windows": ",".join(sorted(common)),
            "total_truth_sites": int(total_truth),
            "total_pred_sites_nonzero": int(total_pred),
            "total_exact_overlap_sites_phased": int(total_phased),
            "total_hamming_errors": int(total_hamming),
            "weighted_hamming_error_rate": weighted_hamming_error_rate,
            "total_switch_errors": int(total_switch),
            "weighted_switch_error_rate": weighted_switch_error_rate,
            "weighted_phased_site_accuracy_pct": weighted_phased_site_accuracy_pct,
            "weighted_truth_correct_pct": weighted_truth_correct_pct,
            "mean_truth_correct_pct": sum(f(r["truth_correct_pct"]) for r in rs) / len(rs),
            "mean_hamming_error_rate": sum(f(r["hamming_error_rate"]) for r in rs) / len(rs),
            "mean_switch_error_rate": sum(f(r["switch_error_rate"]) for r in rs) / len(rs),
            "mean_raw_block_n50_bp": sum(f(r["raw_block_n50_bp"]) for r in rs) / len(rs),
            "mean_median_block_span_bp": sum(f(r["median_block_span_bp"]) for r in rs) / len(rs),
            "mean_max_block_span_bp": sum(f(r["max_block_span_bp"]) for r in rs) / len(rs),
        })

    return out_rows, common


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--per-window-tsv",
        default="results/phase9/phase9a_vs_whatshap.per_window.tsv",
    )
    ap.add_argument(
        "--outdir",
        default="results/phase9",
    )
    ap.add_argument(
        "--prefix",
        default="phase9a_vs_whatshap",
    )
    args = ap.parse_args()

    per_window = Path(args.per_window_tsv)
    outdir = Path(args.outdir)

    rows = read_tsv(per_window)

    print(f"loaded rows: {len(rows)}")
    print(f"input: {per_window}")

    for comparison_set, methods in DEFAULT_METHOD_SETS.items():
        out_rows, common = summarize(rows, methods, comparison_set)
        out_path = outdir / f"{args.prefix}.{comparison_set}.average.tsv"
        write_tsv(out_path, out_rows)

        print()
        print(f"wrote: {out_path}")
        print(f"comparison_set: {comparison_set}")
        print(f"common_windows: {len(common)}")
        if common:
            print(",".join(sorted(common)))

    combined_rows = []
    for comparison_set, methods in DEFAULT_METHOD_SETS.items():
        out_rows, _ = summarize(rows, methods, comparison_set)
        combined_rows.extend(out_rows)

    combined_path = outdir / f"{args.prefix}.all_average_sets.tsv"
    write_tsv(combined_path, combined_rows)

    print()
    print(f"wrote combined table: {combined_path}")


if __name__ == "__main__":
    main()
