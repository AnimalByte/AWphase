#!/usr/bin/env python3
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)


METHODS = [
    {
        "method": "AWPhase_Phase6C_Rust_WMEC",
        "calls": "local_calls.phase6c.tsv",
        "sites": "truth_eval_phase6c.site_comparison.tsv",
        "filled_col": None,
        "prefix": None,
    },
    {
        "method": "AWPhase_Phase8_PBWT",
        "calls": "local_calls.phase8pbwt.tsv",
        "sites": "truth_eval_phase8pbwt.site_comparison.tsv",
        "filled_col": "phase8_pbwt_filled",
        "prefix": "phase8_pbwt",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM",
        "calls": "local_calls.phase8pbwt_hmm.tsv",
        "sites": "truth_eval_phase8pbwt_hmm.site_comparison.tsv",
        "filled_col": "phase8_pbwt_hmm_filled",
        "prefix": "phase8_pbwt_hmm",
    },
    {
        "method": "AWPhase_Phase8_PBWT_V2_exact_prefix",
        "calls": "local_calls.phase8pbwt_v2.tsv",
        "sites": "truth_eval_phase8pbwt_v2.site_comparison.tsv",
        "filled_col": "phase8_pbwt_v2_filled",
        "prefix": "phase8_pbwt_v2",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM_V2_forward_backward",
        "calls": "local_calls.phase8pbwt_hmm_v2.tsv",
        "sites": "truth_eval_phase8pbwt_hmm_v2.site_comparison.tsv",
        "filled_col": "phase8_pbwt_hmm_v2_filled",
        "prefix": "phase8_pbwt_hmm_v2",
    },
    {
        "method": "AWPhase_Phase8_PBWT_V3_bidirectional_prefix",
        "calls": "local_calls.phase8pbwt_v3.tsv",
        "sites": "truth_eval_phase8pbwt_v3.site_comparison.tsv",
        "filled_col": "phase8_pbwt_v3_filled",
        "prefix": "phase8_pbwt_v3",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM_V3_bidirectional_forward_backward",
        "calls": "local_calls.phase8pbwt_hmm_v3.tsv",
        "sites": "truth_eval_phase8pbwt_hmm_v3.site_comparison.tsv",
        "filled_col": "phase8_pbwt_hmm_v3_filled",
        "prefix": "phase8_pbwt_hmm_v3",
    },
]


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


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_json(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with open(path) as handle:
        return json.load(handle)


def load_manifest(path, split):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row = {k: sval(v) for k, v in row.items()}
            if split != "all" and row.get("split") != split:
                continue
            rows.append(row)
    return rows


def source_root(label, source_suffix):
    suffix = source_suffix if source_suffix.startswith("_") else f"_{source_suffix}"
    return Path("results/phase7a_windows") / f"{label}{suffix}"


def is_filled(call_row, method_def):
    col = method_def["filled_col"]
    if not col:
        return False
    return sval(call_row.get(col)).lower() in {"1", "1.0", "true"}


def method_value(call_row, method_def, suffix, default=0.0):
    prefix = method_def["prefix"]
    if prefix and sval(call_row.get(f"{prefix}_{suffix}")) != "":
        return fval(call_row.get(f"{prefix}_{suffix}"), default)
    return fval(call_row.get(suffix), default)


def confidence_bin(confidence):
    x = fval(confidence)
    if x < 0.50:
        return "lt0.50"
    if x < 0.70:
        return "0.50_0.70"
    if x < 0.82:
        return "0.70_0.82"
    if x < 0.90:
        return "0.82_0.90"
    if x < 0.98:
        return "0.90_0.98"
    if x < 1.00:
        return "0.98_1.00"
    return "1.00"


def margin_bin(margin):
    x = fval(margin)
    if x < 1.0:
        return "lt1"
    if x < 5.0:
        return "1_5"
    if x < 25.0:
        return "5_25"
    if x < 100.0:
        return "25_100"
    return "ge100"


def support_bin(value):
    x = fval(value)
    if x < 2:
        return "lt2"
    if x < 4:
        return "2_3"
    if x < 8:
        return "4_7"
    if x < 16:
        return "8_15"
    if x < 32:
        return "16_31"
    return "ge32"


def orientation_by_block(site_rows):
    by_block = defaultdict(list)
    for row in site_rows:
        block = sval(row.get("block_id"))
        pred = ival(row.get("pred_state"), 0)
        truth = ival(row.get("truth_state"), 0)
        matched = sval(row.get("matched_exact_truth")).lower()
        if matched in {"0", "false"} or pred == 0 or truth == 0:
            continue
        if block in {"", "unassigned", "0", ".", "NA", "na"}:
            continue
        by_block[block].append((pred, truth))

    orientation = {}
    counts = {}
    for block, pairs in by_block.items():
        as_is_errors = sum(1 for pred, truth in pairs if pred != truth)
        flipped_errors = sum(1 for pred, truth in pairs if -pred != truth)
        orient = -1 if flipped_errors < as_is_errors else 1
        same = sum(1 for pred, truth in pairs if pred == truth)
        opposite = len(pairs) - same
        orientation[block] = orient
        counts[block] = {
            "orientation_sites": len(pairs),
            "orientation_margin": abs(same - opposite),
            "orientation": orient,
        }
    return orientation, counts


def comparable_sites(site_rows):
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
        if block in {"", "unassigned", "0", ".", "NA", "na"}:
            continue
        orient = orientation.get(block)
        if orient is None:
            continue
        meta = counts.get(block, {})
        out[pos] = {
            "truth_state": truth,
            "oriented_pred_state": pred * orient,
            "correct_after_block_orientation": int(pred * orient == truth),
            "block_orientation": orient,
            "orientation_sites": meta.get("orientation_sites", 0),
            "orientation_margin": meta.get("orientation_margin", 0),
        }
    return out


def calls_by_pos(path):
    out = {}
    for row in read_tsv(path):
        out[ival(row.get("pos"), -1)] = row
    return out


def collect_method_rows(manifest_row, root, method_def):
    calls_path = root / method_def["calls"]
    sites_path = root / method_def["sites"]
    metrics_path = root / method_def["sites"].replace(".site_comparison.tsv", ".metrics.json")
    if not calls_path.exists() or not sites_path.exists():
        return []

    calls = calls_by_pos(calls_path)
    comparable = comparable_sites(read_tsv(sites_path))
    metrics = load_json(metrics_path)
    rows = []

    for pos, detail in comparable.items():
        call = calls.get(pos)
        if not call:
            continue
        state = ival(call.get("local_phase_state") or call.get("phase_state"), 0)
        block = sval(call.get("block_id"))
        if state == 0 or block in {"", "unassigned", "0", ".", "NA", "na"}:
            continue

        filled = is_filled(call, method_def)
        scope = "panel_filled" if filled else "read_backbone"
        confidence = method_value(call, method_def, "confidence", fval(call.get("confidence"), 0.0))
        margin = method_value(call, method_def, "margin", fval(call.get("phase6_margin"), 0.0))
        donors = method_value(call, method_def, "donors", 0.0)
        anchors = method_value(call, method_def, "anchors", fval(call.get("phase6_component_sites"), 0.0))
        read_support = fval(call.get("phase6_site_support"), 0.0)
        read_conflict = fval(call.get("phase6_site_conflict"), 0.0)

        rows.append(
            {
                "split": manifest_row["split"],
                "label": manifest_row["label"],
                "method": method_def["method"],
                "scope": scope,
                "pos": pos,
                "block_id": block,
                "pred_state": state,
                "truth_state": detail["truth_state"],
                "oriented_pred_state": detail["oriented_pred_state"],
                "correct_after_block_orientation": detail["correct_after_block_orientation"],
                "confidence": confidence,
                "margin": margin,
                "donors": donors,
                "anchors": anchors,
                "read_support": read_support,
                "read_conflict": read_conflict,
                "confidence_bin": confidence_bin(confidence),
                "margin_bin": margin_bin(margin),
                "anchors_bin": support_bin(anchors),
                "donors_bin": support_bin(donors),
                "orientation_sites": detail["orientation_sites"],
                "orientation_margin": detail["orientation_margin"],
                "truth_sites_in_window": metrics.get("n_truth_het_sites_in_bed", ""),
            }
        )
    return rows


def aggregate(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["method"], "all_phased", "all", "all")].append(row)
        groups[(row["method"], row["scope"], "all", "all")].append(row)
        for feature in ["confidence_bin", "margin_bin", "anchors_bin", "donors_bin"]:
            groups[(row["method"], "all_phased", feature, row[feature])].append(row)
            groups[(row["method"], row["scope"], feature, row[feature])].append(row)

    out = []
    for (method, scope, feature, bin_name), vals in sorted(groups.items()):
        n = len(vals)
        correct = sum(ival(v["correct_after_block_orientation"]) for v in vals)
        mean_conf = sum(fval(v["confidence"]) for v in vals) / n if n else 0.0
        mean_margin = sum(fval(v["margin"]) for v in vals) / n if n else 0.0
        mean_anchors = sum(fval(v["anchors"]) for v in vals) / n if n else 0.0
        mean_donors = sum(fval(v["donors"]) for v in vals) / n if n else 0.0
        out.append(
            {
                "method": method,
                "scope": scope,
                "feature": feature,
                "bin": bin_name,
                "n_sites": n,
                "n_correct": correct,
                "n_wrong": n - correct,
                "accuracy_pct": 100.0 * correct / n if n else 0.0,
                "error_rate": 1.0 - (correct / n) if n else 0.0,
                "mean_confidence": mean_conf,
                "mean_margin": mean_margin,
                "mean_anchors": mean_anchors,
                "mean_donors": mean_donors,
                "calibration_error": mean_conf - (correct / n) if n else 0.0,
            }
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        default="results/phase8f_manifests/phase8f_windows.window_beds.tsv",
    )
    ap.add_argument("--split", default="train")
    ap.add_argument("--source-suffix", default="_illumina30x_phase7a")
    ap.add_argument("--out-dir", default="results/phase8/confidence_calibration")
    args = ap.parse_args()

    rows = []
    for manifest_row in load_manifest(args.manifest, args.split):
        root = source_root(manifest_row["label"], args.source_suffix)
        if not root.exists():
            continue
        for method_def in METHODS:
            rows.extend(collect_method_rows(manifest_row, root, method_def))

    out_dir = Path(args.out_dir)
    site_fields = [
        "split",
        "label",
        "method",
        "scope",
        "pos",
        "block_id",
        "pred_state",
        "truth_state",
        "oriented_pred_state",
        "correct_after_block_orientation",
        "confidence",
        "margin",
        "donors",
        "anchors",
        "read_support",
        "read_conflict",
        "confidence_bin",
        "margin_bin",
        "anchors_bin",
        "donors_bin",
        "orientation_sites",
        "orientation_margin",
        "truth_sites_in_window",
    ]
    summary_fields = [
        "method",
        "scope",
        "feature",
        "bin",
        "n_sites",
        "n_correct",
        "n_wrong",
        "accuracy_pct",
        "error_rate",
        "mean_confidence",
        "mean_margin",
        "mean_anchors",
        "mean_donors",
        "calibration_error",
    ]

    site_path = out_dir / f"phase8_confidence.{args.split}.sites.tsv"
    summary_path = out_dir / f"phase8_confidence.{args.split}.summary.tsv"
    write_tsv(site_path, rows, site_fields)
    write_tsv(summary_path, aggregate(rows), summary_fields)
    print(f"wrote {site_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
