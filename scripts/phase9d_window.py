#!/usr/bin/env python3
"""Plan and validate current Phase9D benchmark windows.

This wrapper is intentionally data-local: it never downloads or stages external
genomics data. It reads the retained Phase8F/Phase9D window manifest, reports
the current promoted Phase9D default, and checks whether a selected window has
the required local inputs available.
"""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "results/phase8f_manifests/phase8f_windows.window_beds.tsv"
DEFAULT_MODE = "default"
DEFAULT_REF_FASTA = "references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

MODE_THRESHOLDS = {
    "default": "0.80",
    "conservative": "0.90",
    "exploratory": "0.70",
}

REQUIRED_COLUMNS = {
    "split",
    "chrom",
    "start",
    "end",
    "label",
    "bam",
    "panel",
    "json",
    "truth_vcf",
    "truth_bed",
}


@dataclass(frozen=True)
class WindowPlan:
    row: dict[str, str]
    mode: str
    threshold: str
    ref_fasta: str
    output_root: str

    @property
    def label(self) -> str:
        return self.row["label"]

    @property
    def source_window(self) -> str:
        return f"{self.label}_illumina30x_phase7a"

    @property
    def phase7a_root(self) -> str:
        return f"results/phase7a_windows/{self.source_window}"

    @property
    def phase8f_root(self) -> str:
        return f"results/phase8/gated_eval_phase8f/{self.source_window}/t{self.threshold}"

    @property
    def phase8f_calls(self) -> str:
        return f"{self.phase8f_root}/local_calls.phase8f_xgb.tsv"

    @property
    def phase8f_metrics(self) -> str:
        return f"{self.phase8f_root}/truth_eval_phase8f_xgb.metrics.json"

    @property
    def phase9d_candidates(self) -> str:
        return f"{self.phase7a_root}/phase9d_block_bridge_candidates.wide.tsv"

    @property
    def genetic_map(self) -> str:
        return f"data/maps/beagle_grch38/no_chr_in_chrom_field/plink.{self.row['chrom']}.GRCh38.map"

    @property
    def promoted_default_command(self) -> list[str]:
        return [
            "bash",
            "scripts/run_window_awphase_and_whatshap.sh",
            self.row["chrom"],
            self.row["start"],
            self.row["end"],
            self.label,
            self.row["bam"],
            "NONE",
            self.ref_fasta,
        ]

    @property
    def phase7a_command(self) -> list[str]:
        return [
            "bash",
            "scripts/run_phase7a_window.sh",
            self.row["chrom"],
            self.row["start"],
            self.row["end"],
            self.source_window,
            self.row["bam"],
            self.row["panel"],
            self.row["json"],
            self.row["truth_vcf"],
            self.row["truth_bed"],
        ]

    @property
    def phase9d_candidate_command(self) -> list[str]:
        return [
            "python",
            "python/awphase_py/phase9/make_phase9a_block_bridge_candidates_v1.py",
            "--source-window",
            self.source_window,
            "--local-calls-tsv",
            f"{self.phase7a_root}/local_calls.phase6c.tsv",
            "--truth-comparison-tsv",
            f"{self.phase7a_root}/truth_eval_phase6c.site_comparison.tsv",
            "--variant-json",
            f"{self.phase7a_root}/variants.window.json",
            "--panel-bcf",
            self.row["panel"],
            "--genetic-map",
            self.genetic_map,
            "--chrom",
            self.row["chrom"],
            "--start",
            self.row["start"],
            "--end",
            self.row["end"],
            "--max-gap-bp",
            "500000",
            "--max-next-blocks",
            "20",
            "--max-anchors-each-block",
            "16",
            "--min-anchors-each-block",
            "2",
            "--min-orientation-sites",
            "2",
            "--decay-cm",
            "0.05",
            "--top-k",
            "64",
            "--out-tsv",
            self.phase9d_candidates,
        ]


def repo_path(path: str | Path, *, root: Path = ROOT) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return root / candidate


def benchmark_version(path: str | Path) -> str | None:
    text = str(path)
    if "v5.0q" in text or ".v5q." in text:
        return "v5.0q"
    if "v4.2.1" in text:
        return "v4.2.1"
    return None


def load_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"manifest not found: {path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            columns = ", ".join(sorted(missing))
            raise ValueError(f"manifest is missing required columns: {columns}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def select_window(rows: Iterable[dict[str, str]], label: str) -> dict[str, str]:
    matches = [row for row in rows if row.get("label") == label]
    if not matches:
        raise ValueError(f"unknown window label: {label}")
    if len(matches) > 1:
        raise ValueError(f"manifest has duplicate window label: {label}")
    return matches[0]


def sorted_labels(rows: Iterable[dict[str, str]], split: str | None = None) -> list[str]:
    labels = []
    for row in rows:
        if split and row.get("split") != split:
            continue
        labels.append(row["label"])
    return sorted(labels)


def truth_bundle_errors(row: dict[str, str]) -> list[str]:
    vcf_version = benchmark_version(row["truth_vcf"])
    bed_version = benchmark_version(row["truth_bed"])
    if vcf_version and bed_version and vcf_version != bed_version:
        return [
            "truth bundle version mismatch: "
            f"{row['truth_vcf']} is {vcf_version}, {row['truth_bed']} is {bed_version}"
        ]
    return []


def required_input_paths(plan: WindowPlan, *, include_phase9d_inputs: bool = False) -> list[tuple[str, str]]:
    paths = [
        ("illumina_bam", plan.row["bam"]),
        ("panel_bcf", plan.row["panel"]),
        ("variant_json", plan.row["json"]),
        ("truth_vcf", plan.row["truth_vcf"]),
        ("truth_bed", plan.row["truth_bed"]),
        ("reference_fasta", plan.ref_fasta),
    ]
    if include_phase9d_inputs:
        paths.extend(
            [
                ("phase6c_calls", f"{plan.phase7a_root}/local_calls.phase6c.tsv"),
                (
                    "phase6c_truth_comparison",
                    f"{plan.phase7a_root}/truth_eval_phase6c.site_comparison.tsv",
                ),
                ("phase7a_variant_json", f"{plan.phase7a_root}/variants.window.json"),
                ("genetic_map", plan.genetic_map),
            ]
        )
    return paths


def missing_paths(paths: Iterable[tuple[str, str]], *, root: Path = ROOT) -> list[tuple[str, str]]:
    missing = []
    for name, relpath in paths:
        if not repo_path(relpath, root=root).is_file():
            missing.append((name, relpath))
    return missing


def shell_join(command: Iterable[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def render_text_plan(plan: WindowPlan, *, include_phase9d_inputs: bool = False) -> str:
    lines = [
        f"label: {plan.label}",
        f"split: {plan.row['split']}",
        f"region: {plan.row['chrom']}:{plan.row['start']}-{plan.row['end']}",
        "promoted_default: Phase8F_ancestry_weighted_multiwindow",
        f"mode: {plan.mode}",
        f"threshold: {plan.threshold}",
        "",
        "inputs:",
    ]
    for name, relpath in required_input_paths(plan, include_phase9d_inputs=include_phase9d_inputs):
        lines.append(f"  {name}: {relpath}")

    lines.extend(
        [
            "",
            "expected outputs:",
            f"  phase7a_root: {plan.phase7a_root}",
            f"  phase8f_calls: {plan.phase8f_calls}",
            f"  phase8f_metrics: {plan.phase8f_metrics}",
            f"  phase9d_candidates_experimental: {plan.phase9d_candidates}",
            "",
            "commands:",
            "  # broad window runner with WhatsHap Illumina baseline",
            f"  {shell_join(plan.promoted_default_command)}",
            "  # Phase7A prerequisite only",
            f"  {shell_join(plan.phase7a_command)}",
            "  # experimental Phase9D bridge candidate generation",
            f"  PYTHONPATH=python {shell_join(plan.phase9d_candidate_command)}",
        ]
    )
    return "\n".join(lines)


def render_json_plan(plan: WindowPlan, *, include_phase9d_inputs: bool = False) -> str:
    payload = {
        "label": plan.label,
        "split": plan.row["split"],
        "region": {
            "chrom": plan.row["chrom"],
            "start": int(plan.row["start"]),
            "end": int(plan.row["end"]),
        },
        "promoted_default": "Phase8F_ancestry_weighted_multiwindow",
        "mode": plan.mode,
        "threshold": plan.threshold,
        "inputs": dict(required_input_paths(plan, include_phase9d_inputs=include_phase9d_inputs)),
        "outputs": {
            "phase7a_root": plan.phase7a_root,
            "phase8f_calls": plan.phase8f_calls,
            "phase8f_metrics": plan.phase8f_metrics,
            "phase9d_candidates_experimental": plan.phase9d_candidates,
        },
        "commands": {
            "window_runner": plan.promoted_default_command,
            "phase7a_prerequisite": plan.phase7a_command,
            "phase9d_candidate_generation_experimental": plan.phase9d_candidate_command,
        },
    }
    return json.dumps(payload, indent=2, sort_keys=True)


def build_plan(args: argparse.Namespace) -> WindowPlan:
    manifest = repo_path(args.manifest)
    rows = load_manifest(manifest)
    row = select_window(rows, args.label)
    threshold = args.threshold or MODE_THRESHOLDS[args.mode]
    return WindowPlan(
        row=row,
        mode=args.mode,
        threshold=f"{float(threshold):.2f}",
        ref_fasta=args.ref_fasta,
        output_root=args.output_root,
    )


def cmd_list(args: argparse.Namespace) -> int:
    rows = load_manifest(repo_path(args.manifest))
    labels = sorted_labels(rows, args.split)
    if args.format == "json":
        print(json.dumps(labels, indent=2))
    else:
        print("label")
        for label in labels:
            print(label)
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    plan = build_plan(args)
    errors = truth_bundle_errors(plan.row)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    if args.format == "json":
        print(render_json_plan(plan, include_phase9d_inputs=args.include_phase9d_inputs))
    else:
        print(render_text_plan(plan, include_phase9d_inputs=args.include_phase9d_inputs))
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    plan = build_plan(args)
    errors = truth_bundle_errors(plan.row)
    missing = missing_paths(
        required_input_paths(plan, include_phase9d_inputs=args.include_phase9d_inputs)
    )

    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    for name, relpath in missing:
        print(f"MISSING {name}: {relpath}", file=sys.stderr)

    if errors or missing:
        return 1

    print(f"{plan.label}: inputs present")
    return 0


def add_window_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("label", help="window label from the manifest, such as chr22_20_25mb")
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST.relative_to(ROOT)),
        help="Phase8F/Phase9D window manifest TSV",
    )
    parser.add_argument(
        "--mode",
        choices=sorted(MODE_THRESHOLDS),
        default=DEFAULT_MODE,
        help="promoted Phase9D mode",
    )
    parser.add_argument(
        "--threshold",
        help="override the mode threshold; default comes from --mode",
    )
    parser.add_argument(
        "--ref-fasta",
        default=DEFAULT_REF_FASTA,
        help="locally staged GRCh38 reference FASTA path",
    )
    parser.add_argument(
        "--output-root",
        default="results/phase9d_windows",
        help="reserved output root for future canonical Phase9D runs",
    )
    parser.add_argument(
        "--include-phase9d-inputs",
        action="store_true",
        help="include experimental Phase9D bridge-generation prerequisites",
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="list manifest window labels")
    list_parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST.relative_to(ROOT)),
        help="Phase8F/Phase9D window manifest TSV",
    )
    list_parser.add_argument("--split", help="optional split filter, such as train")
    list_parser.add_argument("--format", choices=["text", "json"], default="text")
    list_parser.set_defaults(func=cmd_list)

    plan_parser = subparsers.add_parser("plan", help="print a no-download run plan")
    add_window_args(plan_parser)
    plan_parser.add_argument("--format", choices=["text", "json"], default="text")
    plan_parser.set_defaults(func=cmd_plan)

    check_parser = subparsers.add_parser("check", help="validate local inputs for a window")
    add_window_args(check_parser)
    check_parser.set_defaults(func=cmd_check)

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
