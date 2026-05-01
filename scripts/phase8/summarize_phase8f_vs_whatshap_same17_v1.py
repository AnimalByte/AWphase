#!/usr/bin/env python3
import csv
import json
from pathlib import Path
from collections import defaultdict

MANIFEST = Path("results/phase8f_manifests/phase8f_windows.window_beds.tsv")
P8_SUMMARY = Path("results/phase8/gated_eval_phase8f/phase8f_gate_eval.summary.tsv")

OUT_PER = Path("results/phase8/phase8f_vs_whatshap.same17.per_window.tsv")
OUT_AVG = Path("results/phase8/phase8f_vs_whatshap.same17.average_by_available.tsv")
OUT_COMMON_ALL = Path("results/phase8/phase8f_vs_whatshap.same17.common_all_methods.average.tsv")
OUT_COMMON_ILL = Path("results/phase8/phase8f_vs_whatshap.same17.common_illumina.average.tsv")
OUT_COMMON_HIFI = Path("results/phase8/phase8f_vs_whatshap.same17.common_hifi.average.tsv")

METRICS = [
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
]

def clean_row(r):
    return {k: (v.replace("\r", "") if isinstance(v, str) else v) for k, v in r.items()}

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def add_json(rows, label, method, path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return
    m = load_json(p)
    row = {"label": label, "method": method, "path": str(p)}
    for k in METRICS:
        row[k] = m.get(k, "")
    rows.append(row)

def add_phase8f(rows, labels, threshold):
    p8 = {}
    if not P8_SUMMARY.exists():
        return
    for r in csv.DictReader(open(P8_SUMMARY), delimiter="\t"):
        if r["threshold"] == threshold:
            p8[r["source_window"]] = r

    for label in labels:
        source_window = f"{label}_illumina30x_phase7a"
        if source_window not in p8:
            continue
        r = p8[source_window]
        row = {
            "label": label,
            "method": f"AWPhase_Phase8F_t{threshold}",
            "path": f"results/phase8/gated_eval_phase8f/{source_window}/t{threshold}/truth_eval_phase8f_xgb.metrics.json",
        }
        for k in METRICS:
            row[k] = r.get(k, "")
        # Phase8F summary table does not have these two denominators/truth-sites.
        # Pull them from the metrics JSON if available.
        mp = Path(row["path"])
        if mp.exists():
            m = load_json(mp)
            for k in METRICS:
                row[k] = m.get(k, row.get(k, ""))
        rows.append(row)

def f(x):
    try:
        return float(x)
    except Exception:
        return 0.0

def summarize(rows, out, method_filter=None, label_filter=None, comparison_set="available"):
    rs = rows
    if method_filter is not None:
        rs = [r for r in rs if r["method"] in method_filter]
    if label_filter is not None:
        rs = [r for r in rs if r["label"] in label_filter]

    by = defaultdict(list)
    for r in rs:
        by[r["method"]].append(r)

    fields = [
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

    out_rows = []
    for method, mr in sorted(by.items()):
        if not mr:
            continue

        total_truth = sum(f(r["n_truth_het_sites_in_bed"]) for r in mr)
        total_phased = sum(f(r["n_exact_overlap_sites_phased"]) for r in mr)
        total_hamm_den = sum(f(r["hamming_denominator"]) for r in mr)
        total_switch_den = sum(f(r["switch_denominator"]) for r in mr)
        total_hamm = sum(f(r["hamming_errors"]) for r in mr)
        total_switch = sum(f(r["switch_errors"]) for r in mr)

        out_rows.append({
            "comparison_set": comparison_set,
            "method": method,
            "n_windows": len(mr),
            "windows": ",".join(r["label"] for r in mr),
            "total_truth_sites": int(total_truth),
            "total_pred_sites_nonzero": int(sum(f(r["n_pred_sites_nonzero"]) for r in mr)),
            "total_exact_overlap_sites_phased": int(total_phased),
            "total_hamming_errors": int(total_hamm),
            "weighted_hamming_error_rate": total_hamm / total_hamm_den if total_hamm_den else 0,
            "total_switch_errors": int(total_switch),
            "weighted_switch_error_rate": total_switch / total_switch_den if total_switch_den else 0,
            "weighted_phased_site_accuracy_pct": 100.0 * (1.0 - total_hamm / total_hamm_den) if total_hamm_den else 0,
            "weighted_truth_correct_pct": 100.0 * (total_phased - total_hamm) / total_truth if total_truth else 0,
            "mean_truth_correct_pct": sum(f(r["truth_correct_pct"]) for r in mr) / len(mr),
            "mean_hamming_error_rate": sum(f(r["hamming_error_rate"]) for r in mr) / len(mr),
            "mean_switch_error_rate": sum(f(r["switch_error_rate"]) for r in mr) / len(mr),
            "mean_raw_block_n50_bp": sum(f(r["raw_block_n50_bp"]) for r in mr) / len(mr),
            "mean_median_block_span_bp": sum(f(r["median_block_span_bp"]) for r in mr) / len(mr),
        })

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print("wrote", out)

labels = []
for r in csv.DictReader(open(MANIFEST), delimiter="\t"):
    r = clean_row(r)
    if r["split"] == "train":
        labels.append(r["label"])

rows = []

for label in labels:
    root = Path(f"results/phase7a_windows/{label}_illumina30x_phase7a")
    add_json(rows, label, "AWPhase_Phase6C_read_backbone", root / "truth_eval_phase6c.metrics.json")
    add_json(rows, label, "AWPhase_Phase7A_hand_panel_fill", root / "truth_eval_phase7a.metrics.json")
    add_json(rows, label, "WhatsHap_Illumina30x", Path(f"baselines/whatshap/{label}_whatshap_illumina30x/truth_eval.metrics.json"))
    add_json(rows, label, "WhatsHap_GIAB_CCS_HiFi30x", Path(f"baselines/whatshap/{label}_whatshap_hifi30x/truth_eval.metrics.json"))

add_phase8f(rows, labels, "0.50")
add_phase8f(rows, labels, "0.80")
add_phase8f(rows, labels, "0.90")

fields = ["label", "method"] + METRICS + ["path"]
OUT_PER.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_PER, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("wrote", OUT_PER)

summarize(rows, OUT_AVG, comparison_set="all_available_per_method")

# Common Illumina set: methods requiring Illumina baseline.
ill_methods = [
    "AWPhase_Phase6C_read_backbone",
    "AWPhase_Phase7A_hand_panel_fill",
    "AWPhase_Phase8F_t0.50",
    "AWPhase_Phase8F_t0.80",
    "AWPhase_Phase8F_t0.90",
    "WhatsHap_Illumina30x",
]
labels_by_method = defaultdict(set)
for r in rows:
    labels_by_method[r["method"]].add(r["label"])

common_ill = set(labels)
for m in ill_methods:
    common_ill &= labels_by_method[m]

summarize(rows, OUT_COMMON_ILL, method_filter=ill_methods, label_filter=common_ill, comparison_set="common_with_illumina_whatshap")

hifi_methods = [
    "AWPhase_Phase6C_read_backbone",
    "AWPhase_Phase7A_hand_panel_fill",
    "AWPhase_Phase8F_t0.50",
    "AWPhase_Phase8F_t0.80",
    "AWPhase_Phase8F_t0.90",
    "WhatsHap_GIAB_CCS_HiFi30x",
]

common_hifi = set(labels)
for m in hifi_methods:
    common_hifi &= labels_by_method[m]

summarize(rows, OUT_COMMON_HIFI, method_filter=hifi_methods, label_filter=common_hifi, comparison_set="common_with_hifi_whatshap")

all_methods = sorted(set(ill_methods + hifi_methods))
common_all = set(labels)
for m in all_methods:
    common_all &= labels_by_method[m]

summarize(rows, OUT_COMMON_ALL, method_filter=all_methods, label_filter=common_all, comparison_set="common_all_methods")

print("common_ill", sorted(common_ill))
print("common_hifi", sorted(common_hifi))
print("common_all", sorted(common_all))
