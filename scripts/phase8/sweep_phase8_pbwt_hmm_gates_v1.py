#!/usr/bin/env python3
import argparse
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
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


@dataclass(frozen=True)
class GateConfig:
    name: str
    min_confidence: float
    min_anchors: int
    min_donors: int
    pbwt_min_margin: float
    hmm_min_margin: float
    min_best_vs_second_margin: float
    allow_hmm_override: bool


@dataclass(frozen=True)
class Candidate:
    pos: int
    block_id: str
    pred_state: int
    support_plus: float
    support_minus: float
    margin: float
    confidence: float
    donors: int
    anchors: int
    method: str
    best_vs_second_margin: float = 0.0
    n_block_candidates: int = 0


DEFAULT_CONFIGS = [
    GateConfig("current_like", 0.82, 2, 12, 1.0, 1.0, 0.25, True),
    GateConfig("pbwt_only_current", 0.82, 2, 12, 1.0, 1.0, 0.25, False),
    GateConfig("anchor4_donor32", 0.82, 4, 32, 1.0, 1.0, 0.25, True),
    GateConfig("conf90_anchor4_margin10", 0.90, 4, 32, 10.0, 1.0, 1.0, True),
    GateConfig("conf98_anchor4_margin25", 0.98, 4, 32, 25.0, 1.0, 5.0, True),
    GateConfig("pbwt_high_margin", 0.90, 4, 64, 100.0, 1.0, 25.0, True),
    GateConfig("hmm_soft_override", 0.90, 4, 32, 10.0, 0.75, 0.25, True),
    GateConfig("hmm_soft_anchor8", 0.90, 8, 32, 10.0, 0.75, 0.25, True),
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
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_with_headers(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [dict(r) for r in reader]
        return list(reader.fieldnames or []), rows


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


def parse_candidate(row):
    return Candidate(
        pos=ival(row.get("pos"), -1),
        block_id=sval(row.get("block_id")),
        pred_state=ival(row.get("pred_state"), 0),
        support_plus=fval(row.get("support_plus"), 0.0),
        support_minus=fval(row.get("support_minus"), 0.0),
        margin=fval(row.get("margin"), 0.0),
        confidence=fval(row.get("confidence"), 0.0),
        donors=ival(row.get("donors"), 0),
        anchors=ival(row.get("anchors"), 0),
        method=sval(row.get("method")),
    )


def rank_best(candidates, config):
    if not candidates:
        return None
    ranked = sorted(
        candidates,
        key=lambda c: (
            -c.margin,
            -c.confidence,
            -c.donors,
            -c.anchors,
            c.block_id,
            c.pred_state,
        ),
    )
    best = ranked[0]
    second_margin = ranked[1].margin if len(ranked) > 1 else 0.0
    best = Candidate(
        pos=best.pos,
        block_id=best.block_id,
        pred_state=best.pred_state,
        support_plus=best.support_plus,
        support_minus=best.support_minus,
        margin=best.margin,
        confidence=best.confidence,
        donors=best.donors,
        anchors=best.anchors,
        method=best.method,
        best_vs_second_margin=best.margin - second_margin,
        n_block_candidates=len(ranked),
    )
    min_margin = config.hmm_min_margin if best.method == "pbwt_hmm" else config.pbwt_min_margin
    if best.donors < config.min_donors:
        return None
    if best.anchors < config.min_anchors:
        return None
    if best.confidence < config.min_confidence:
        return None
    if best.margin < min_margin:
        return None
    if best.best_vs_second_margin < config.min_best_vs_second_margin:
        return None
    return best


def choose_candidates(candidate_rows, config):
    by_pos_method = {}
    for row in candidate_rows:
        cand = parse_candidate(row)
        if cand.pos < 0 or cand.pred_state == 0 or not cand.block_id:
            continue
        if cand.method not in {"pbwt", "pbwt_hmm"}:
            continue
        by_pos_method.setdefault((cand.pos, cand.method), []).append(cand)

    chosen = {}
    for (pos, method), candidates in by_pos_method.items():
        if method != "pbwt":
            continue
        best = rank_best(candidates, config)
        if best is not None:
            chosen[pos] = best

    if config.allow_hmm_override:
        for (pos, method), candidates in by_pos_method.items():
            if method != "pbwt_hmm":
                continue
            best = rank_best(candidates, config)
            if best is not None:
                chosen[pos] = best

    return chosen


def ensure_headers(headers):
    out = list(headers)
    for field in [
        "phase8_pbwt_hmm_filled",
        "phase8_pbwt_hmm_confidence",
        "phase8_pbwt_hmm_margin",
        "phase8_pbwt_hmm_support_plus",
        "phase8_pbwt_hmm_support_minus",
        "phase8_pbwt_hmm_donors",
        "phase8_pbwt_hmm_anchors",
        "phase8_pbwt_hmm_method",
        "phase8_pbwt_hmm_best_vs_second_margin",
        "phase8_pbwt_hmm_n_block_candidates",
    ]:
        if field not in out:
            out.append(field)
    return out


def apply_gates(backbone_tsv, candidates_tsv, out_tsv, config):
    headers, rows = read_tsv_with_headers(backbone_tsv)
    candidate_rows = read_tsv(candidates_tsv)
    chosen = choose_candidates(candidate_rows, config)
    headers = ensure_headers(headers)

    out_rows = []
    for row in rows:
        out = dict(row)
        pos = ival(out.get("pos"), -1)
        pred = chosen.get(pos)
        if pred is not None:
            if "local_phase_state" in out:
                out["local_phase_state"] = str(pred.pred_state)
            if "phase_state" in out:
                out["phase_state"] = str(pred.pred_state)
            out["block_id"] = pred.block_id
            out["phase8_pbwt_hmm_filled"] = "1"
            out["phase8_pbwt_hmm_confidence"] = f"{pred.confidence:.6f}"
            out["phase8_pbwt_hmm_margin"] = f"{pred.margin:.6f}"
            out["phase8_pbwt_hmm_support_plus"] = f"{pred.support_plus:.6f}"
            out["phase8_pbwt_hmm_support_minus"] = f"{pred.support_minus:.6f}"
            out["phase8_pbwt_hmm_donors"] = str(pred.donors)
            out["phase8_pbwt_hmm_anchors"] = str(pred.anchors)
            out["phase8_pbwt_hmm_method"] = pred.method
            out["phase8_pbwt_hmm_best_vs_second_margin"] = f"{pred.best_vs_second_margin:.6f}"
            out["phase8_pbwt_hmm_n_block_candidates"] = str(pred.n_block_candidates)
        else:
            out.setdefault("phase8_pbwt_hmm_filled", "0")
            if out.get("phase8_pbwt_hmm_filled") != "1":
                out["phase8_pbwt_hmm_filled"] = "0"
                for field in headers:
                    if field.startswith("phase8_pbwt_hmm_") and field != "phase8_pbwt_hmm_filled":
                        out.setdefault(field, "")
        out_rows.append(out)

    write_tsv(out_tsv, out_rows, headers)
    return len(chosen)


def run_eval_task(task):
    manifest_row, config_dict, out_dir = task
    config = GateConfig(**config_dict)
    label = manifest_row["label"]
    split = manifest_row["split"]
    root = result_root(label)
    run_dir = Path(out_dir) / split / config.name / label
    run_dir.mkdir(parents=True, exist_ok=True)

    backbone_tsv = root / "local_calls.phase6c.tsv"
    candidates_tsv = root / "phase8_pbwt_hmm" / "candidates.pbwt_hmm.tsv"
    variant_json = root / "variants.window.json"
    out_calls = run_dir / "local_calls.phase8pbwt_hmm.gated.tsv"
    out_prefix = run_dir / "truth_eval"

    for path in [backbone_tsv, candidates_tsv, variant_json, Path(manifest_row["truth_vcf"]), Path(manifest_row["truth_bed"])]:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing required file for {label}: {path}")

    accepted_sites = apply_gates(backbone_tsv, candidates_tsv, out_calls, config)
    env = dict(os.environ)
    env["PYTHONPATH"] = "python"
    subprocess.run(
        [
            sys.executable,
            "python/awphase_py/evaluate_phase_truth.py",
            "--pred-tsv",
            str(out_calls),
            "--variant-json",
            str(variant_json),
            "--truth-vcf",
            manifest_row["truth_vcf"],
            "--truth-bed",
            manifest_row["truth_bed"],
            "--sample",
            "HG002",
            "--phase-column",
            "local_phase_state",
            "--out-prefix",
            str(out_prefix),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    with open(str(out_prefix) + ".metrics.json") as fh:
        metrics = json.load(fh)
    out = {
        "split": split,
        "label": label,
        "config": config.name,
        "accepted_sites_by_gate": accepted_sites,
        "calls_tsv": str(out_calls),
        "metrics_json": str(out_prefix) + ".metrics.json",
    }
    out.update(asdict(config))
    for key in METRICS:
        out[key] = metrics.get(key, "")
    return out


def evaluate_configs(manifest_rows, configs, out_dir, jobs):
    tasks = [
        (row, asdict(config), str(out_dir))
        for row in manifest_rows
        for config in configs
    ]
    results = []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futures = [ex.submit(run_eval_task, task) for task in tasks]
        for i, fut in enumerate(as_completed(futures), start=1):
            results.append(fut.result())
            if i % 10 == 0 or i == len(tasks):
                print(f"completed {i}/{len(tasks)} gate evaluations", file=sys.stderr)
    results.sort(key=lambda r: (r["config"], r["split"], r["label"]))
    return results


def summarize(rows):
    by_config = {}
    for row in rows:
        by_config.setdefault(row["config"], []).append(row)

    summary = []
    for config_name, config_rows in sorted(by_config.items()):
        truth = sum(fval(r.get("n_truth_het_sites_in_bed")) for r in config_rows)
        phased = sum(fval(r.get("n_exact_overlap_sites_phased")) for r in config_rows)
        hden = sum(fval(r.get("hamming_denominator")) for r in config_rows)
        herr = sum(fval(r.get("hamming_errors")) for r in config_rows)
        sden = sum(fval(r.get("switch_denominator")) for r in config_rows)
        serr = sum(fval(r.get("switch_errors")) for r in config_rows)
        correct = phased - herr
        template = config_rows[0]
        weighted_hamming = herr / hden if hden else 0.0
        weighted_switch = serr / sden if sden else 0.0
        weighted_truth_correct = 100.0 * correct / truth if truth else 0.0
        summary.append(
            {
                "config": config_name,
                "n_windows": len(config_rows),
                "windows": ",".join(r["label"] for r in config_rows),
                "min_confidence": template["min_confidence"],
                "min_anchors": template["min_anchors"],
                "min_donors": template["min_donors"],
                "pbwt_min_margin": template["pbwt_min_margin"],
                "hmm_min_margin": template["hmm_min_margin"],
                "min_best_vs_second_margin": template["min_best_vs_second_margin"],
                "allow_hmm_override": template["allow_hmm_override"],
                "total_truth_sites": int(truth),
                "total_pred_sites_nonzero": int(sum(fval(r.get("n_pred_sites_nonzero")) for r in config_rows)),
                "total_exact_overlap_sites_phased": int(phased),
                "total_accepted_sites_by_gate": int(sum(fval(r.get("accepted_sites_by_gate")) for r in config_rows)),
                "total_hamming_errors": int(herr),
                "weighted_hamming_error_rate": weighted_hamming,
                "total_switch_errors": int(serr),
                "weighted_switch_error_rate": weighted_switch,
                "weighted_phased_site_accuracy_pct": 100.0 * (1.0 - weighted_hamming) if hden else 0.0,
                "weighted_truth_correct_pct": weighted_truth_correct,
                "mean_truth_correct_pct": sum(fval(r.get("truth_correct_pct")) for r in config_rows) / len(config_rows),
                "mean_hamming_error_rate": sum(fval(r.get("hamming_error_rate")) for r in config_rows) / len(config_rows),
                "mean_switch_error_rate": sum(fval(r.get("switch_error_rate")) for r in config_rows) / len(config_rows),
                "mean_raw_block_n50_bp": sum(fval(r.get("raw_block_n50_bp")) for r in config_rows) / len(config_rows),
            }
        )

    current = next((r for r in summary if r["config"] == "current_like"), None)
    for row in summary:
        if current:
            row["delta_truth_correct_vs_current"] = (
                row["weighted_truth_correct_pct"] - current["weighted_truth_correct_pct"]
            )
            row["delta_hamming_vs_current"] = (
                row["weighted_hamming_error_rate"] - current["weighted_hamming_error_rate"]
            )
            row["delta_switch_vs_current"] = (
                row["weighted_switch_error_rate"] - current["weighted_switch_error_rate"]
            )
        else:
            row["delta_truth_correct_vs_current"] = ""
            row["delta_hamming_vs_current"] = ""
            row["delta_switch_vs_current"] = ""
    summary.sort(
        key=lambda r: (
            -fval(r["weighted_truth_correct_pct"]),
            fval(r["weighted_hamming_error_rate"]),
            fval(r["weighted_switch_error_rate"]),
        )
    )
    return summary


def select_train_config(summary_rows):
    current = next((r for r in summary_rows if r["config"] == "current_like"), None)
    if current is None:
        return summary_rows[0]["config"] if summary_rows else ""
    max_hamming = fval(current["weighted_hamming_error_rate"]) * 1.05
    max_switch = fval(current["weighted_switch_error_rate"]) * 1.10
    eligible = [
        r
        for r in summary_rows
        if fval(r["weighted_hamming_error_rate"]) <= max_hamming
        and fval(r["weighted_switch_error_rate"]) <= max_switch
    ]
    if not eligible:
        return "current_like"
    eligible.sort(
        key=lambda r: (
            -fval(r["weighted_truth_correct_pct"]),
            fval(r["weighted_hamming_error_rate"]),
            fval(r["weighted_switch_error_rate"]),
        )
    )
    return eligible[0]["config"]


def config_by_name(name):
    by_name = {c.name: c for c in DEFAULT_CONFIGS}
    return by_name[name]


def output_fields(per_window=False):
    config_fields = [
        "config",
        "min_confidence",
        "min_anchors",
        "min_donors",
        "pbwt_min_margin",
        "hmm_min_margin",
        "min_best_vs_second_margin",
        "allow_hmm_override",
    ]
    if per_window:
        return [
            "split",
            "label",
        ] + config_fields + [
            "accepted_sites_by_gate",
        ] + METRICS + [
            "calls_tsv",
            "metrics_json",
        ]
    return config_fields + [
        "n_windows",
        "windows",
        "total_truth_sites",
        "total_pred_sites_nonzero",
        "total_exact_overlap_sites_phased",
        "total_accepted_sites_by_gate",
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
        "delta_truth_correct_vs_current",
        "delta_hamming_vs_current",
        "delta_switch_vs_current",
        "selected_by_train",
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        default="results/phase8f_manifests/phase8f_windows.window_beds.tsv",
    )
    ap.add_argument("--out-dir", default="results/phase8/pbwt_hmm_gate_sweep")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    train_rows = load_manifest(args.manifest, "train")
    holdout_rows = load_manifest(args.manifest, "holdout_reserved")

    print(f"running train sweep: {len(DEFAULT_CONFIGS)} configs x {len(train_rows)} windows", file=sys.stderr)
    train_results = evaluate_configs(train_rows, DEFAULT_CONFIGS, out_dir, args.jobs)
    train_summary = summarize(train_results)
    selected = select_train_config(train_summary)
    for row in train_summary:
        row["selected_by_train"] = int(row["config"] == selected)

    write_tsv(out_dir / "train.per_window.tsv", train_results, output_fields(per_window=True))
    write_tsv(out_dir / "train.summary.tsv", train_summary, output_fields(per_window=False))

    holdout_configs = [config_by_name("current_like")]
    if selected != "current_like":
        holdout_configs.append(config_by_name(selected))

    print(
        f"selected {selected} from train; running holdout for {', '.join(c.name for c in holdout_configs)}",
        file=sys.stderr,
    )
    holdout_results = evaluate_configs(holdout_rows, holdout_configs, out_dir, args.jobs)
    holdout_summary = summarize(holdout_results)
    for row in holdout_summary:
        row["selected_by_train"] = int(row["config"] == selected)

    write_tsv(out_dir / "holdout.per_window.tsv", holdout_results, output_fields(per_window=True))
    write_tsv(out_dir / "holdout.summary.tsv", holdout_summary, output_fields(per_window=False))

    print(f"wrote {out_dir / 'train.summary.tsv'}")
    print(f"wrote {out_dir / 'holdout.summary.tsv'}")
    print(f"selected_by_train={selected}")


if __name__ == "__main__":
    main()
