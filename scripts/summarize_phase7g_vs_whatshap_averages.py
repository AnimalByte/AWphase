#!/usr/bin/env python3
import json
import csv
from pathlib import Path
from statistics import mean

LABELS = [
    "chr1_15_20mb",
    "chr20_25_30mb",
    "chr20_30_35mb",
    "chr20_45_50mb",
    "chr22_15_20mb",
    "chr22_20_25mb",
]

OUTDIR = Path("results/summary")
OUTDIR.mkdir(parents=True, exist_ok=True)

def load_json(path):
    p = Path(path)
    if not p.exists():
        return None
    with open(p) as fh:
        return json.load(fh)

def add(rows, label, method, path):
    m = load_json(path)
    if m is None:
        return

    rows.append({
        "label": label,
        "method": method,
        "path": str(path),
        "n_truth_het_sites_in_bed": m.get("n_truth_het_sites_in_bed", 0) or 0,
        "n_pred_sites_nonzero": m.get("n_pred_sites_nonzero", 0) or 0,
        "n_exact_overlap_sites_phased": m.get("n_exact_overlap_sites_phased", 0) or 0,
        "hamming_errors": m.get("hamming_errors", 0) or 0,
        "hamming_denominator": m.get("hamming_denominator", m.get("n_exact_overlap_sites_phased", 0)) or 0,
        "hamming_error_rate": m.get("hamming_error_rate", 0) or 0,
        "switch_errors": m.get("switch_errors", 0) or 0,
        "switch_denominator": m.get("switch_denominator", 0) or 0,
        "switch_error_rate": m.get("switch_error_rate", 0) or 0,
        "phased_site_accuracy_pct": m.get("phased_site_accuracy_pct", 0) or 0,
        "truth_correct_pct": m.get("truth_correct_pct", 0) or 0,
        "raw_block_n50_bp": m.get("raw_block_n50_bp", 0) or 0,
        "median_block_span_bp": m.get("median_block_span_bp", 0) or 0,
    })

rows = []

for label in LABELS:
    root = Path(f"results/phase7a_windows/{label}_illumina30x_phase7a")

    # AWPhase components
    add(rows, label, "AWPhase_Phase6C_read_backbone", root / "truth_eval_phase6c.metrics.json")
    add(rows, label, "AWPhase_Phase7A_hand_panel_fill", root / "truth_eval_phase7a.metrics.json")

    add(
        rows,
        label,
        "AWPhase_Phase7G_donor_default_t0.50",
        Path(f"results/phase7e_xgboost_acceptance/gated_eval_phase7g/{label}_illumina30x_phase7a/t0.50/truth_eval_phase7g_xgb.metrics.json")
    )
    add(
        rows,
        label,
        "AWPhase_Phase7G_donor_precision_t0.80",
        Path(f"results/phase7e_xgboost_acceptance/gated_eval_phase7g/{label}_illumina30x_phase7a/t0.80/truth_eval_phase7g_xgb.metrics.json")
    )

    # WhatsHap Illumina 30x
    add(
        rows,
        label,
        "WhatsHap_Illumina30x",
        Path(f"baselines/whatshap/{label}_whatshap_illumina30x/truth_eval.metrics.json")
    )

    # WhatsHap HiFi/PacBio baseline; discover possible naming variants.
    for name in [
        f"baselines/whatshap/{label}_whatshap_hifi30x/truth_eval.metrics.json",
        f"baselines/whatshap/{label}_whatshap_hifi/truth_eval.metrics.json",
        f"baselines/whatshap/{label}_whatshap_pacbio30x/truth_eval.metrics.json",
        f"baselines/whatshap/{label}_whatshap_pacbio/truth_eval.metrics.json",
    ]:
        if Path(name).exists():
            add(rows, label, "WhatsHap_HiFi_or_PacBio30x", Path(name))
            break

# Write per-window table.
per_window = OUTDIR / "phase7g_vs_whatshap.per_window.tsv"
fields = [
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
    "path",
]

with open(per_window, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)

def summarize(subset, name, out_path):
    methods = sorted(set(r["method"] for r in subset))
    out_rows = []

    for method in methods:
        rs = [r for r in subset if r["method"] == method]
        if not rs:
            continue

        hden = sum(float(r["hamming_denominator"]) for r in rs)
        herr = sum(float(r["hamming_errors"]) for r in rs)
        sden = sum(float(r["switch_denominator"]) for r in rs)
        serr = sum(float(r["switch_errors"]) for r in rs)

        truth = sum(float(r["n_truth_het_sites_in_bed"]) for r in rs)
        phased = sum(float(r["n_exact_overlap_sites_phased"]) for r in rs)
        correct = phased - herr

        weighted_hamming = herr / hden if hden else 0.0
        weighted_switch = serr / sden if sden else 0.0
        weighted_site_acc = 100.0 * (1.0 - weighted_hamming) if hden else 0.0
        weighted_truth_correct = 100.0 * correct / truth if truth else 0.0

        out_rows.append({
            "comparison_set": name,
            "method": method,
            "n_windows": len(rs),
            "windows": ",".join(sorted(r["label"] for r in rs)),
            "total_truth_sites": int(truth),
            "total_pred_sites_nonzero": int(sum(float(r["n_pred_sites_nonzero"]) for r in rs)),
            "total_exact_overlap_sites_phased": int(phased),
            "total_hamming_errors": int(herr),
            "weighted_hamming_error_rate": weighted_hamming,
            "total_switch_errors": int(serr),
            "weighted_switch_error_rate": weighted_switch,
            "weighted_phased_site_accuracy_pct": weighted_site_acc,
            "weighted_truth_correct_pct": weighted_truth_correct,
            "mean_truth_correct_pct": mean(float(r["truth_correct_pct"]) for r in rs),
            "mean_hamming_error_rate": mean(float(r["hamming_error_rate"]) for r in rs),
            "mean_switch_error_rate": mean(float(r["switch_error_rate"]) for r in rs),
            "mean_raw_block_n50_bp": mean(float(r["raw_block_n50_bp"]) for r in rs),
            "mean_median_block_span_bp": mean(float(r["median_block_span_bp"]) for r in rs),
        })

    out_fields = [
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
    ]

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for r in sorted(out_rows, key=lambda x: (x["comparison_set"], x["method"])):
            w.writerow(r)

    return out_rows

# 1. All available windows per method.
summarize(rows, "all_available_windows_per_method", OUTDIR / "phase7g_vs_whatshap.average_all_available.tsv")

# 2. Fair comparison against WhatsHap Illumina: only labels where Illumina WhatsHap exists.
illumina_labels = sorted(set(
    r["label"] for r in rows if r["method"] == "WhatsHap_Illumina30x"
))
illumina_subset = [
    r for r in rows
    if r["label"] in illumina_labels
    and r["method"] in {
        "AWPhase_Phase6C_read_backbone",
        "AWPhase_Phase7A_hand_panel_fill",
        "AWPhase_Phase7G_donor_default_t0.50",
        "AWPhase_Phase7G_donor_precision_t0.80",
        "WhatsHap_Illumina30x",
    }
]
summarize(illumina_subset, "fair_common_windows_with_whatshap_illumina30x", OUTDIR / "phase7g_vs_whatshap.average_common_illumina.tsv")

# 3. Fair comparison against HiFi/PacBio: only labels where HiFi/PacBio WhatsHap exists.
hifi_labels = sorted(set(
    r["label"] for r in rows if r["method"] == "WhatsHap_HiFi_or_PacBio30x"
))
hifi_subset = [
    r for r in rows
    if r["label"] in hifi_labels
    and r["method"] in {
        "AWPhase_Phase6C_read_backbone",
        "AWPhase_Phase7A_hand_panel_fill",
        "AWPhase_Phase7G_donor_default_t0.50",
        "AWPhase_Phase7G_donor_precision_t0.80",
        "WhatsHap_HiFi_or_PacBio30x",
    }
]
summarize(hifi_subset, "fair_common_windows_with_whatshap_hifi_or_pacbio30x", OUTDIR / "phase7g_vs_whatshap.average_common_hifi.tsv")

print("wrote", per_window)
print("wrote", OUTDIR / "phase7g_vs_whatshap.average_all_available.tsv")
print("wrote", OUTDIR / "phase7g_vs_whatshap.average_common_illumina.tsv")
print("wrote", OUTDIR / "phase7g_vs_whatshap.average_common_hifi.tsv")
print()
print("illumina_common_windows:", illumina_labels)
print("hifi_common_windows:", hifi_labels)
