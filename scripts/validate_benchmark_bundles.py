#!/usr/bin/env python3
import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SHELL_PATHS = [
    ROOT / "scripts",
]
TSV_PATHS = [
    ROOT / "results/phase8f_manifests/phase8f_windows.tsv",
    ROOT / "results/phase8f_manifests/phase8f_windows.window_beds.tsv",
]


def display_path(path):
    try:
        return path.relative_to(ROOT)
    except ValueError:
        return path


def benchmark_version(path):
    text = str(path)
    if "v5.0q" in text or ".v5q." in text:
        return "v5.0q"
    if "v4.2.1" in text:
        return "v4.2.1"
    return None


def extract_shell_value(line):
    match = re.search(r"=\$?\{?[A-Z_]*:?-?([^}\"']+)", line)
    if match:
        return match.group(1)
    match = re.search(r"=[\"']?([^\"'\n]+)", line)
    if match:
        return match.group(1)
    return line


def validate_shell_literal_pairs(errors):
    for root in SHELL_PATHS:
        for path in sorted(root.rglob("*.sh")):
            lines = path.read_text().splitlines()
            pending_vcf = None
            for lineno, line in enumerate(lines, start=1):
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue

                if re.search(r"\bTRUTH_VCF=", stripped):
                    value = extract_shell_value(stripped)
                    version = benchmark_version(value)
                    if version:
                        pending_vcf = (lineno, value, version)
                    continue

                if re.search(r"\bTRUTH_BED=", stripped):
                    value = extract_shell_value(stripped)
                    bed_version = benchmark_version(value)
                    if not bed_version or pending_vcf is None:
                        continue
                    vcf_lineno, vcf_value, vcf_version = pending_vcf
                    if lineno - vcf_lineno > 8:
                        continue
                    if vcf_version != bed_version:
                        rel = display_path(path)
                        errors.append(
                            f"{rel}:{vcf_lineno}-{lineno}: mixed benchmark bundle "
                            f"{vcf_version} VCF with {bed_version} BED "
                            f"({vcf_value} vs {value})"
                        )


def validate_manifest_pairs(errors):
    for path in TSV_PATHS:
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if "truth_vcf" not in reader.fieldnames or "truth_bed" not in reader.fieldnames:
                continue
            for lineno, row in enumerate(reader, start=2):
                vcf = row.get("truth_vcf", "")
                bed = row.get("truth_bed", "")
                vcf_version = benchmark_version(vcf)
                bed_version = benchmark_version(bed)
                if vcf_version and bed_version and vcf_version != bed_version:
                    rel = display_path(path)
                    label = row.get("label", "unknown")
                    errors.append(
                        f"{rel}:{lineno}: {label} mixes {vcf_version} VCF with "
                        f"{bed_version} BED ({vcf} vs {bed})"
                    )


def main():
    errors = []
    validate_shell_literal_pairs(errors)
    validate_manifest_pairs(errors)

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("benchmark bundle versions ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
