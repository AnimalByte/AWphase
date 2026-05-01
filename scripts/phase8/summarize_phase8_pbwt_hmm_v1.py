#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

METRICS = [
    "n_truth_het_sites_in_bed",
    "n_pred_sites_nonzero",
    "n_exact_overlap_sites_phased",
    "pct_truth_sites_phased",
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
]


def f(value):
    try:
        return float(value)
    except Exception:
        return 0.0


def median(values):
    values = sorted(values)
    if not values:
        return 0.0
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2.0


def population_stddev(values):
    if not values:
        return 0.0
    mean = sum(values) / len(values)
    variance = sum((v - mean) ** 2 for v in values) / len(values)
    return variance ** 0.5


def load_manifest(path, split):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row = {k: (v or "").strip() for k, v in row.items()}
            if split != "all" and row["split"] != split:
                continue
            rows.append(row)
    return rows


def load_runtime_seconds(path):
    path = Path(path)
    out = {}
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            method = (row.get("method") or "").strip()
            if not method:
                continue
            out[method] = out.get(method, 0.0) + f(row.get("runtime_seconds", 0))
    return out


def source_root(label, source_suffix):
    suffix = source_suffix if source_suffix.startswith("_") else f"_{source_suffix}"
    return Path(f"results/phase7a_windows/{label}{suffix}")


def add_json(rows, split, label, method, path, runtime_seconds=""):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return
    with open(path) as fh:
        metrics = json.load(fh)
    row = {"split": split, "label": label, "method": method, "path": str(path)}
    for key in METRICS:
        row[key] = metrics.get(key, "")
    truth = f(metrics.get("n_truth_het_sites_in_bed", 0))
    phased = f(metrics.get("n_exact_overlap_sites_phased", 0))
    row["phasing_coverage_pct"] = 100.0 * phased / truth if truth else 0.0
    row["runtime_seconds"] = runtime_seconds
    rows.append(row)


def summarize(rows):
    by_method = defaultdict(list)
    for row in rows:
        by_method[row["method"]].append(row)

    out = []
    for method, method_rows in sorted(by_method.items()):
        truth = sum(f(r["n_truth_het_sites_in_bed"]) for r in method_rows)
        phased = sum(f(r["n_exact_overlap_sites_phased"]) for r in method_rows)
        hden = sum(f(r["hamming_denominator"]) for r in method_rows)
        herr = sum(f(r["hamming_errors"]) for r in method_rows)
        sden = sum(f(r["switch_denominator"]) for r in method_rows)
        serr = sum(f(r["switch_errors"]) for r in method_rows)
        correct = phased - herr
        raw_n50_values = [f(r["raw_block_n50_bp"]) for r in method_rows]
        runtime_values = [
            f(r["runtime_seconds"])
            for r in method_rows
            if str(r.get("runtime_seconds", "")).strip() != ""
        ]

        out.append(
            {
                "method": method,
                "n_windows": len(method_rows),
                "windows": ",".join(r["label"] for r in method_rows),
                "total_truth_sites": int(truth),
                "total_pred_sites_nonzero": int(
                    sum(f(r["n_pred_sites_nonzero"]) for r in method_rows)
                ),
                "total_exact_overlap_sites_phased": int(phased),
                "weighted_phasing_coverage_pct": 100.0 * phased / truth if truth else 0.0,
                "total_hamming_errors": int(herr),
                "weighted_hamming_error_rate": herr / hden if hden else 0.0,
                "total_switch_errors": int(serr),
                "weighted_switch_error_rate": serr / sden if sden else 0.0,
                "weighted_phased_site_accuracy_pct": 100.0 * (1.0 - herr / hden)
                if hden
                else 0.0,
                "weighted_truth_correct_pct": 100.0 * correct / truth if truth else 0.0,
                "mean_truth_correct_pct": sum(
                    f(r["truth_correct_pct"]) for r in method_rows
                )
                / len(method_rows),
                "mean_hamming_error_rate": sum(
                    f(r["hamming_error_rate"]) for r in method_rows
                )
                / len(method_rows),
                "mean_switch_error_rate": sum(
                    f(r["switch_error_rate"]) for r in method_rows
                )
                / len(method_rows),
                "mean_raw_block_n50_bp": sum(
                    f(r["raw_block_n50_bp"]) for r in method_rows
                )
                / len(method_rows),
                "median_raw_block_n50_bp": median(raw_n50_values),
                "stddev_raw_block_n50_bp": population_stddev(raw_n50_values),
                "total_runtime_seconds": sum(runtime_values) if runtime_values else 0.0,
                "mean_runtime_seconds": sum(runtime_values) / len(runtime_values)
                if runtime_values
                else 0.0,
                "median_runtime_seconds": median(runtime_values),
                "stddev_runtime_seconds": population_stddev(runtime_values),
            }
        )
    return out


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        default="results/phase8f_manifests/phase8f_windows.window_beds.tsv",
    )
    ap.add_argument("--split", default="train")
    ap.add_argument(
        "--out-dir",
        default="results/phase8/pbwt_hmm_manifest",
    )
    ap.add_argument("--source-suffix", default="_illumina30x_phase7a")
    args = ap.parse_args()

    manifest_rows = load_manifest(args.manifest, args.split)
    rows = []
    for row in manifest_rows:
        split = row["split"]
        label = row["label"]
        root = source_root(label, args.source_suffix)
        runtimes = {}
        runtimes.update(load_runtime_seconds(root / "run_timing.tsv"))
        runtimes.update(load_runtime_seconds(root / "phase8_timing.tsv"))
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase6C_Rust_WMEC",
            root / "truth_eval_phase6c.metrics.json",
            runtimes.get("AWPhase_Phase6C_Rust_WMEC", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase7A_current",
            root / "truth_eval_phase7a.metrics.json",
            runtimes.get("AWPhase_Phase7A_current", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT",
            root / "truth_eval_phase8pbwt.metrics.json",
            runtimes.get("AWPhase_Phase8_PBWT", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT_HMM",
            root / "truth_eval_phase8pbwt_hmm.metrics.json",
            runtimes.get("AWPhase_Phase8_PBWT_HMM", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT_V2_exact_prefix",
            root / "truth_eval_phase8pbwt_v2.metrics.json",
            runtimes.get("AWPhase_Phase8_PBWT_V2_exact_prefix", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT_HMM_V2_forward_backward",
            root / "truth_eval_phase8pbwt_hmm_v2.metrics.json",
            runtimes.get("AWPhase_Phase8_PBWT_HMM_V2_forward_backward", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT_V3_bidirectional_prefix",
            root / "truth_eval_phase8pbwt_v3.metrics.json",
            runtimes.get("AWPhase_Phase8_PBWT_V3_bidirectional_prefix", ""),
        )
        add_json(
            rows,
            split,
            label,
            "AWPhase_Phase8_PBWT_HMM_V3_bidirectional_forward_backward",
            root / "truth_eval_phase8pbwt_hmm_v3.metrics.json",
            runtimes.get(
                "AWPhase_Phase8_PBWT_HMM_V3_bidirectional_forward_backward",
                "",
            ),
        )
        illumina_runtime = load_runtime_seconds(
            Path(f"baselines/whatshap/{label}_whatshap_illumina30x/run_timing.tsv")
        ).get("WhatsHap_Illumina30x", "")
        add_json(
            rows,
            split,
            label,
            "WhatsHap_Illumina30x",
            Path(f"baselines/whatshap/{label}_whatshap_illumina30x/truth_eval.metrics.json"),
            illumina_runtime,
        )
        pacbio_runtime = load_runtime_seconds(
            Path(f"baselines/whatshap/{label}_whatshap_pacbio/run_timing.tsv")
        ).get("WhatsHap_PacBio_HiFi30x", "")
        add_json(
            rows,
            split,
            label,
            "WhatsHap_PacBio_HiFi30x",
            Path(f"baselines/whatshap/{label}_whatshap_pacbio/truth_eval.metrics.json"),
            pacbio_runtime,
        )

    out_dir = Path(args.out_dir)
    per_window = out_dir / f"phase8_pbwt_hmm.{args.split}.per_window.tsv"
    average = out_dir / f"phase8_pbwt_hmm.{args.split}.average.tsv"

    write_tsv(
        per_window,
        rows,
        ["split", "label", "method"]
        + METRICS
        + ["phasing_coverage_pct", "runtime_seconds", "path"],
    )
    avg_rows = summarize(rows)
    avg_fields = [
        "method",
        "n_windows",
        "windows",
        "total_truth_sites",
        "total_pred_sites_nonzero",
        "total_exact_overlap_sites_phased",
        "weighted_phasing_coverage_pct",
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
        "median_raw_block_n50_bp",
        "stddev_raw_block_n50_bp",
        "total_runtime_seconds",
        "mean_runtime_seconds",
        "median_runtime_seconds",
        "stddev_runtime_seconds",
    ]
    write_tsv(average, avg_rows, avg_fields)

    print(f"wrote {per_window}")
    print(f"wrote {average}")


if __name__ == "__main__":
    main()
