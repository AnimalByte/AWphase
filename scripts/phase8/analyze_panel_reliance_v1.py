#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


csv.field_size_limit(sys.maxsize)

METHODS = [
    {
        "method": "AWPhase_Phase6C_Rust_WMEC",
        "calls": "local_calls.phase6c.tsv",
        "filled_col": None,
        "summary": "solve.phase6c.summary.json",
        "metrics": "truth_eval_phase6c.metrics.json",
        "sites": "truth_eval_phase6c.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase7A_current",
        "calls": "local_calls.phase7a.tsv",
        "filled_col": "phase7a_panel_filled",
        "summary": "fill.phase7a.summary.json",
        "metrics": "truth_eval_phase7a.metrics.json",
        "sites": "truth_eval_phase7a.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT",
        "calls": "local_calls.phase8pbwt.tsv",
        "filled_col": "phase8_pbwt_filled",
        "summary": "phase8_pbwt/summary.pbwt.json",
        "metrics": "truth_eval_phase8pbwt.metrics.json",
        "sites": "truth_eval_phase8pbwt.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM",
        "calls": "local_calls.phase8pbwt_hmm.tsv",
        "filled_col": "phase8_pbwt_hmm_filled",
        "summary": "phase8_pbwt_hmm/summary.pbwt_hmm.json",
        "metrics": "truth_eval_phase8pbwt_hmm.metrics.json",
        "sites": "truth_eval_phase8pbwt_hmm.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT_V2_exact_prefix",
        "calls": "local_calls.phase8pbwt_v2.tsv",
        "filled_col": "phase8_pbwt_v2_filled",
        "summary": "phase8_pbwt_v2/summary.pbwt_v2.json",
        "metrics": "truth_eval_phase8pbwt_v2.metrics.json",
        "sites": "truth_eval_phase8pbwt_v2.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM_V2_forward_backward",
        "calls": "local_calls.phase8pbwt_hmm_v2.tsv",
        "filled_col": "phase8_pbwt_hmm_v2_filled",
        "summary": "phase8_pbwt_hmm_v2/summary.pbwt_hmm_v2.json",
        "metrics": "truth_eval_phase8pbwt_hmm_v2.metrics.json",
        "sites": "truth_eval_phase8pbwt_hmm_v2.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT_V3_bidirectional_prefix",
        "calls": "local_calls.phase8pbwt_v3.tsv",
        "filled_col": "phase8_pbwt_v3_filled",
        "summary": "phase8_pbwt_v3/summary.pbwt_v3.json",
        "metrics": "truth_eval_phase8pbwt_v3.metrics.json",
        "sites": "truth_eval_phase8pbwt_v3.site_comparison.tsv",
    },
    {
        "method": "AWPhase_Phase8_PBWT_HMM_V3_bidirectional_forward_backward",
        "calls": "local_calls.phase8pbwt_hmm_v3.tsv",
        "filled_col": "phase8_pbwt_hmm_v3_filled",
        "summary": "phase8_pbwt_hmm_v3/summary.pbwt_hmm_v3.json",
        "metrics": "truth_eval_phase8pbwt_hmm_v3.metrics.json",
        "sites": "truth_eval_phase8pbwt_hmm_v3.site_comparison.tsv",
    },
]


NUMERIC_FIELDS = [
    "truth_sites",
    "variant_sites",
    "read_backbone_phased_sites",
    "method_phased_sites",
    "method_filled_sites",
    "method_panel_added_sites",
    "method_panel_added_with_read_evidence_sites",
    "method_panel_added_without_read_evidence_sites",
    "method_preserved_read_sites",
    "method_changed_read_sites",
    "method_lost_read_sites",
    "panel_exact_variant_sites",
    "panel_missing_variant_sites",
    "panel_allele_mismatch_variant_sites",
    "panel_nonbiallelic_variant_sites",
    "panel_needed_sites",
    "panel_found_sites",
    "panel_missing_needed_sites",
    "panel_allele_mismatch_needed_records",
    "panel_nonbiallelic_needed_records",
    "candidate_unphased_positions",
    "accepted_sites",
    "panel_added_panel_exact_sites",
    "panel_added_panel_missing_sites",
    "panel_added_panel_allele_mismatch_sites",
    "panel_added_panel_nonbiallelic_sites",
]


def f(value):
    try:
        return float(value)
    except Exception:
        return 0.0


def load_manifest(path, split):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row = {k: (v or "").strip() for k, v in row.items()}
            if split != "all" and row["split"] != split:
                continue
            rows.append(row)
    return rows


def load_json(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with open(path) as handle:
        return json.load(handle)


def load_call_rows(path):
    out = {}
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pos = int(float(row["pos"]))
            out[pos] = row
    return out


def load_site_rows(path):
    out = {}
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("matched_exact_truth", "1")).strip() in {"0", "false", "False"}:
                continue
            pos = int(float(row["pos"]))
            out[pos] = row
    return out


def is_phased(row):
    if not row:
        return False
    state = int(float(row.get("local_phase_state") or row.get("phase_state") or 0))
    block_id = (row.get("block_id") or "").strip()
    return state != 0 and block_id not in {"", "unassigned", "0", ".", "NA"}


def is_eval_phased(row):
    if not row:
        return False
    state = int(float(row.get("pred_state") or row.get("local_phase_state") or 0))
    row_block = (row.get("block_id") or "").strip()
    return state != 0 and row_block not in {"", "unassigned", "0", ".", "NA"}


def phase_state(row):
    if not row:
        return 0
    return int(float(row.get("local_phase_state") or row.get("phase_state") or 0))


def eval_phase_state(row):
    if not row:
        return 0
    return int(float(row.get("pred_state") or row.get("local_phase_state") or 0))


def block_id(row):
    return (row.get("block_id") or "").strip()


def read_evidence(row):
    if not row:
        return 0.0
    return f(row.get("phase6_site_support")) + f(row.get("phase6_site_conflict"))


def filled(row, col):
    if not row or not col:
        return False
    return str(row.get(col, "")).strip() in {"1", "1.0", "true", "True"}


def benchmark_variants(site_rows):
    out = {}
    for pos, row in site_rows.items():
        out[pos] = (
            str(row.get("truth_ref") or row.get("ref") or "").upper(),
            str(row.get("truth_alt") or row.get("alt") or "").upper(),
        )
    return out


def load_panel_status(panel_bcf, chrom, start, end, variants):
    statuses = {pos: "panel_missing" for pos in variants}
    cmd = [
        "bcftools",
        "query",
        "-r",
        f"{chrom}:{start}-{end}",
        "-f",
        "%POS\t%REF\t%ALT\n",
        panel_bcf,
    ]
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        pos_s, ref, alt = line.split("\t")[:3]
        pos = int(pos_s)
        if pos not in variants:
            continue
        if "," in alt:
            if statuses[pos] != "panel_exact":
                statuses[pos] = "panel_nonbiallelic"
            continue
        ref = ref.upper()
        alt = alt.upper()
        expected_ref, expected_alt = variants[pos]
        if ref == expected_ref and alt == expected_alt:
            statuses[pos] = "panel_exact"
        elif statuses[pos] == "panel_missing":
            statuses[pos] = "panel_allele_mismatch"
    return statuses


def method_paths(root, method_def):
    return {
        "calls": root / method_def["calls"],
        "summary": root / method_def["summary"],
        "metrics": root / method_def["metrics"],
        "sites": root / method_def["sites"],
    }


def source_root(label, source_suffix):
    suffix = source_suffix if source_suffix.startswith("_") else f"_{source_suffix}"
    return Path(f"results/phase7a_windows/{label}{suffix}")


def summarize_method(
    split,
    manifest_row,
    method_def,
    base_call_rows,
    method_call_rows,
    base_site_rows,
    method_site_rows,
    variants,
    panel_status,
    root,
):
    paths = method_paths(root, method_def)
    metrics = load_json(paths["metrics"])
    summary = load_json(paths["summary"])

    eval_positions = set(variants)
    read_phased = {
        pos for pos, row in base_site_rows.items() if pos in eval_positions and is_eval_phased(row)
    }
    method_phased = {
        pos
        for pos, row in method_site_rows.items()
        if pos in eval_positions and is_eval_phased(row)
    }
    method_filled = {
        pos
        for pos, row in method_call_rows.items()
        if pos in eval_positions and filled(row, method_def["filled_col"])
    }
    panel_added = method_phased - read_phased
    panel_added_with_read = {
        pos for pos in panel_added if read_evidence(base_call_rows.get(pos)) > 0.0
    }
    panel_added_without_read = panel_added - panel_added_with_read

    preserved_read = {
        pos
        for pos in read_phased & method_phased
        if eval_phase_state(base_site_rows[pos]) == eval_phase_state(method_site_rows[pos])
        and block_id(base_site_rows[pos]) == block_id(method_site_rows[pos])
    }
    changed_read = (read_phased & method_phased) - preserved_read
    lost_read = read_phased - method_phased

    status_counts = defaultdict(int)
    for status in panel_status.values():
        status_counts[status] += 1

    added_status_counts = defaultdict(int)
    for pos in panel_added:
        added_status_counts[panel_status.get(pos, "panel_missing")] += 1

    panel_needed = (
        summary.get("panel_positions_needed")
        or summary.get("needed_panel_positions")
        or 0
    )
    panel_found = summary.get("panel_positions_found", 0)
    panel_missing = summary.get("panel_positions_missing", 0)

    return {
        "split": split,
        "label": manifest_row["label"],
        "method": method_def["method"],
        "truth_sites": int(f(metrics.get("n_truth_het_sites_in_bed"))),
        "variant_sites": len(variants),
        "read_backbone_phased_sites": len(read_phased),
        "method_phased_sites": len(method_phased),
        "method_filled_sites": len(method_filled),
        "method_panel_added_sites": len(panel_added),
        "method_panel_added_with_read_evidence_sites": len(panel_added_with_read),
        "method_panel_added_without_read_evidence_sites": len(panel_added_without_read),
        "method_preserved_read_sites": len(preserved_read),
        "method_changed_read_sites": len(changed_read),
        "method_lost_read_sites": len(lost_read),
        "panel_added_pct_truth_sites": 100.0
        * len(panel_added)
        / f(metrics.get("n_truth_het_sites_in_bed"))
        if f(metrics.get("n_truth_het_sites_in_bed"))
        else 0.0,
        "panel_added_pct_method_phased_sites": 100.0
        * len(panel_added)
        / len(method_phased)
        if method_phased
        else 0.0,
        "panel_exact_variant_sites": status_counts["panel_exact"],
        "panel_missing_variant_sites": status_counts["panel_missing"],
        "panel_allele_mismatch_variant_sites": status_counts["panel_allele_mismatch"],
        "panel_nonbiallelic_variant_sites": status_counts["panel_nonbiallelic"],
        "panel_exact_pct_variant_sites": 100.0
        * status_counts["panel_exact"]
        / len(variants)
        if variants
        else 0.0,
        "panel_needed_sites": int(f(panel_needed)),
        "panel_found_sites": int(f(panel_found)),
        "panel_missing_needed_sites": int(f(panel_missing)),
        "panel_found_pct_needed_sites": 100.0 * f(panel_found) / f(panel_needed)
        if f(panel_needed)
        else 0.0,
        "panel_allele_mismatch_needed_records": int(
            f(summary.get("panel_allele_mismatch_records"))
        ),
        "panel_nonbiallelic_needed_records": int(
            f(summary.get("panel_nonbiallelic_records"))
        ),
        "candidate_unphased_positions": int(
            f(summary.get("candidate_unphased_positions"))
        ),
        "accepted_sites": int(
            f(summary.get("accepted_sites") or summary.get("filled_sites"))
        ),
        "panel_added_panel_exact_sites": added_status_counts["panel_exact"],
        "panel_added_panel_missing_sites": added_status_counts["panel_missing"],
        "panel_added_panel_allele_mismatch_sites": added_status_counts[
            "panel_allele_mismatch"
        ],
        "panel_added_panel_nonbiallelic_sites": added_status_counts[
            "panel_nonbiallelic"
        ],
    }


def aggregate(rows):
    by_method = defaultdict(list)
    for row in rows:
        by_method[row["method"]].append(row)

    out = []
    for method, method_rows in sorted(by_method.items()):
        agg = {"method": method, "n_windows": len(method_rows)}
        agg["windows"] = ",".join(row["label"] for row in method_rows)
        for field in NUMERIC_FIELDS:
            agg[field] = sum(f(row.get(field)) for row in method_rows)

        truth = agg["truth_sites"]
        method_phased = agg["method_phased_sites"]
        variants = agg["variant_sites"]
        panel_needed = agg["panel_needed_sites"]
        agg["weighted_panel_added_pct_truth_sites"] = (
            100.0 * agg["method_panel_added_sites"] / truth if truth else 0.0
        )
        agg["weighted_panel_added_pct_method_phased_sites"] = (
            100.0 * agg["method_panel_added_sites"] / method_phased
            if method_phased
            else 0.0
        )
        agg["weighted_panel_exact_pct_variant_sites"] = (
            100.0 * agg["panel_exact_variant_sites"] / variants if variants else 0.0
        )
        agg["weighted_panel_found_pct_needed_sites"] = (
            100.0 * agg["panel_found_sites"] / panel_needed if panel_needed else 0.0
        )
        out.append(agg)
    return out


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--split", default="train")
    ap.add_argument("--out-dir", default="results/phase8/panel_reliance_manifest")
    ap.add_argument("--source-suffix", default="_illumina30x_phase7a")
    args = ap.parse_args()

    manifest_rows = load_manifest(args.manifest, args.split)
    rows = []

    for manifest_row in manifest_rows:
        chrom = manifest_row["chrom"]
        start = int(manifest_row["start"])
        end = int(manifest_row["end"])
        label = manifest_row["label"]
        root = source_root(label, args.source_suffix)
        base_call_rows = load_call_rows(root / "local_calls.phase6c.tsv")
        base_site_rows = load_site_rows(root / "truth_eval_phase6c.site_comparison.tsv")
        if not base_call_rows or not base_site_rows:
            continue

        variants = benchmark_variants(base_site_rows)
        panel_status = load_panel_status(
            manifest_row["panel"],
            chrom,
            start,
            end,
            variants,
        )

        for method_def in METHODS:
            paths = method_paths(root, method_def)
            method_call_rows = load_call_rows(paths["calls"])
            method_site_rows = load_site_rows(paths["sites"])
            if not method_call_rows or not method_site_rows:
                continue
            rows.append(
                summarize_method(
                    manifest_row["split"],
                    manifest_row,
                    method_def,
                    base_call_rows,
                    method_call_rows,
                    base_site_rows,
                    method_site_rows,
                    variants,
                    panel_status,
                    root,
                )
            )

    out_dir = Path(args.out_dir)
    per_window = out_dir / f"panel_reliance.{args.split}.per_window.tsv"
    average = out_dir / f"panel_reliance.{args.split}.average.tsv"

    per_fields = [
        "split",
        "label",
        "method",
        *NUMERIC_FIELDS,
        "panel_added_pct_truth_sites",
        "panel_added_pct_method_phased_sites",
        "panel_exact_pct_variant_sites",
        "panel_found_pct_needed_sites",
    ]
    avg_fields = [
        "method",
        "n_windows",
        "windows",
        *NUMERIC_FIELDS,
        "weighted_panel_added_pct_truth_sites",
        "weighted_panel_added_pct_method_phased_sites",
        "weighted_panel_exact_pct_variant_sites",
        "weighted_panel_found_pct_needed_sites",
    ]
    write_tsv(per_window, rows, per_fields)
    write_tsv(average, aggregate(rows), avg_fields)

    print(f"wrote {per_window}")
    print(f"wrote {average}")


if __name__ == "__main__":
    main()
