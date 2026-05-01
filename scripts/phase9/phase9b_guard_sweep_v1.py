#!/usr/bin/env python3
import csv
import json
import subprocess
from pathlib import Path
from collections import defaultdict

MANIFEST = Path("results/phase8f_manifests/phase8f_windows.window_beds.tsv")
PRED = Path("results/phase9/phase9a_bridge_xgb.noleak.oof_predictions.tsv")
OUTROOT = Path("results/phase9/phase9b_guard_sweep")
BASE_PHASE8F = Path("results/phase8/gated_eval_phase8f")

CONFIGS = [
    # name, min_prob, min_conf, min_margin, min_eff_donors, min_usable, max_gap_bp, min_total_support
    ("p995_c990_m25_e2_u2_g100k_s25", 0.995, 0.990, 25, 2, 2, 100000, 25),
    ("p995_c995_m25_e2_u2_g100k_s25", 0.995, 0.995, 25, 2, 2, 100000, 25),
    ("p995_c995_m50_e3_u2_g75k_s50", 0.995, 0.995, 50, 3, 2, 75000, 50),
    ("p995_c995_m100_e4_u2_g50k_s100", 0.995, 0.995, 100, 4, 2, 50000, 100),
    ("p998_c995_m50_e3_u2_g75k_s50", 0.998, 0.995, 50, 3, 2, 75000, 50),
    ("p998_c997_m100_e4_u2_g50k_s100", 0.998, 0.997, 100, 4, 2, 50000, 100),
    ("p999_c997_m100_e4_u2_g50k_s100", 0.999, 0.997, 100, 4, 2, 50000, 100),
    ("p999_c999_m200_e6_u3_g25k_s200", 0.999, 0.999, 200, 6, 3, 25000, 200),
]

def clean(v):
    return v.replace("\r", "") if isinstance(v, str) else v

def fval(x, default=0.0):
    try:
        if x is None or str(x).strip() == "":
            return default
        return float(x)
    except Exception:
        return default

def load_manifest():
    rows = []
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k: clean(v) for k, v in r.items()}
            if r["split"] == "train":
                rows.append(r)
    return rows

def filter_predictions(config):
    name, min_prob, min_conf, min_margin, min_eff, min_usable, max_gap, min_support = config
    out = OUTROOT / "filtered_predictions" / f"{name}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    total = 0

    with open(PRED) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        with open(out, "w", newline="") as oh:
            w = csv.DictWriter(oh, fieldnames=fields, delimiter="\t")
            w.writeheader()

            for r in reader:
                total += 1

                if fval(r.get("xgb_oof_probability")) < min_prob:
                    continue
                if fval(r.get("bridge_confidence")) < min_conf:
                    continue
                if fval(r.get("bridge_margin")) < min_margin:
                    continue
                if fval(r.get("bridge_effective_donors")) < min_eff:
                    continue
                if fval(r.get("bridge_a_usable_anchors")) < min_usable:
                    continue
                if fval(r.get("bridge_b_usable_anchors")) < min_usable:
                    continue
                if fval(r.get("gap_bp")) > max_gap:
                    continue
                if fval(r.get("bridge_total_support")) < min_support:
                    continue

                w.writerow(r)
                kept += 1

    return out, kept, total

def run_cmd(cmd, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    return p.returncode

def main():
    manifest = load_manifest()
    OUTROOT.mkdir(parents=True, exist_ok=True)

    average_rows = []
    per_window_all = []

    for config in CONFIGS:
        name, min_prob, *_ = config
        filtered_pred, kept, total = filter_predictions(config)

        summary = OUTROOT / name / "phase9b_guard_eval.summary.tsv"
        summary.parent.mkdir(parents=True, exist_ok=True)

        with open(summary, "w", newline="") as sh:
            fields = [
                "source_window", "config", "n_pred_sites_nonzero",
                "hamming_errors", "hamming_error_rate",
                "switch_errors", "switch_error_rate",
                "phased_site_accuracy_pct", "truth_correct_pct",
                "raw_block_n50_bp", "median_block_span_bp", "max_block_span_bp",
                "accepted_prediction_rows"
            ]
            wsum = csv.DictWriter(sh, fieldnames=fields, delimiter="\t")
            wsum.writeheader()

            for r in manifest:
                label = r["label"]
                source_window = f"{label}_illumina30x_phase7a"
                root = Path(f"results/phase7a_windows/{source_window}")
                base_calls = BASE_PHASE8F / source_window / "t0.80" / "local_calls.phase8f_xgb.tsv"

                if not base_calls.exists():
                    continue

                outdir = OUTROOT / name / source_window
                outdir.mkdir(parents=True, exist_ok=True)

                local_out = outdir / "local_calls.phase9b_bridge.tsv"
                apply_log = outdir / "apply_bridge.log"
                eval_log = outdir / "eval.log"
                eval_prefix = outdir / "truth_eval_phase9b_bridge"

                cmd = [
                    "python", "python/awphase_py/phase9/apply_phase9a_block_bridges_v1.py",
                    "--local-calls-tsv", str(base_calls),
                    "--bridge-predictions-tsv", str(filtered_pred),
                    "--source-window", source_window,
                    "--threshold", str(min_prob),
                    "--out-tsv", str(local_out),
                ]
                rc = run_cmd(cmd, apply_log)
                if rc != 0:
                    continue

                cmd = [
                    "python", "python/awphase_py/evaluate_phase_truth.py",
                    "--pred-tsv", str(local_out),
                    "--variant-json", str(root / "variants.window.json"),
                    "--truth-vcf", r["truth_vcf"],
                    "--truth-bed", r["truth_bed"],
                    "--phase-column", "local_phase_state",
                    "--out-prefix", str(eval_prefix),
                ]
                rc = run_cmd(cmd, eval_log)
                if rc != 0:
                    continue

                m = json.load(open(str(eval_prefix) + ".metrics.json"))
                row = {
                    "source_window": source_window,
                    "config": name,
                    "n_pred_sites_nonzero": m.get("n_pred_sites_nonzero"),
                    "hamming_errors": m.get("hamming_errors"),
                    "hamming_error_rate": m.get("hamming_error_rate"),
                    "switch_errors": m.get("switch_errors"),
                    "switch_error_rate": m.get("switch_error_rate"),
                    "phased_site_accuracy_pct": m.get("phased_site_accuracy_pct"),
                    "truth_correct_pct": m.get("truth_correct_pct"),
                    "raw_block_n50_bp": m.get("raw_block_n50_bp"),
                    "median_block_span_bp": m.get("median_block_span_bp"),
                    "max_block_span_bp": m.get("max_block_span_bp"),
                    "accepted_prediction_rows": kept,
                }
                wsum.writerow(row)
                per_window_all.append(row)

        rows = list(csv.DictReader(open(summary), delimiter="\t"))

        def ff(x):
            try:
                return float(x)
            except Exception:
                return 0.0

        if rows:
            average_rows.append({
                "config": name,
                "n_windows": len(rows),
                "accepted_prediction_rows": kept,
                "total_pred_sites_nonzero": int(sum(ff(r["n_pred_sites_nonzero"]) for r in rows)),
                "total_hamming_errors": int(sum(ff(r["hamming_errors"]) for r in rows)),
                "mean_hamming_error_rate": sum(ff(r["hamming_error_rate"]) for r in rows) / len(rows),
                "total_switch_errors": int(sum(ff(r["switch_errors"]) for r in rows)),
                "mean_switch_error_rate": sum(ff(r["switch_error_rate"]) for r in rows) / len(rows),
                "mean_phased_site_accuracy_pct": sum(ff(r["phased_site_accuracy_pct"]) for r in rows) / len(rows),
                "mean_truth_correct_pct": sum(ff(r["truth_correct_pct"]) for r in rows) / len(rows),
                "mean_raw_block_n50_bp": sum(ff(r["raw_block_n50_bp"]) for r in rows) / len(rows),
                "mean_median_block_span_bp": sum(ff(r["median_block_span_bp"]) for r in rows) / len(rows),
                "mean_max_block_span_bp": sum(ff(r["max_block_span_bp"]) for r in rows) / len(rows),
            })

    avg_out = OUTROOT / "phase9b_guard_sweep.average.tsv"
    avg_fields = [
        "config", "n_windows", "accepted_prediction_rows",
        "total_pred_sites_nonzero", "total_hamming_errors", "mean_hamming_error_rate",
        "total_switch_errors", "mean_switch_error_rate",
        "mean_phased_site_accuracy_pct", "mean_truth_correct_pct",
        "mean_raw_block_n50_bp", "mean_median_block_span_bp", "mean_max_block_span_bp",
    ]

    with open(avg_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=avg_fields, delimiter="\t")
        w.writeheader()
        w.writerows(average_rows)

    print("wrote", avg_out)

if __name__ == "__main__":
    main()
