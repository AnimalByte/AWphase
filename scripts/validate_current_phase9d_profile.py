#!/usr/bin/env python3
import csv
import json
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "configs/current_phase9d.profile.yaml"


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def require_path(relpath, errors):
    path = ROOT / relpath
    if not path.exists():
        errors.append(f"missing path: {relpath}")
    return path


def main():
    errors = []

    if not PROFILE.exists():
        return fail(f"missing profile: {PROFILE.relative_to(ROOT)}")

    with PROFILE.open() as fh:
        profile = yaml.safe_load(fh)

    if not isinstance(profile, dict):
        return fail("profile did not parse as a mapping")

    status = profile.get("status", {})
    if status.get("current_project_target") != "Phase9D":
        errors.append("status.current_project_target must be Phase9D")
    if status.get("promoted_default") != "Phase8F_ancestry_weighted_multiwindow":
        errors.append("status.promoted_default must be Phase8F_ancestry_weighted_multiwindow")

    expected_thresholds = {
        "default": 0.80,
        "conservative": 0.90,
        "exploratory": 0.70,
    }
    modes = profile.get("modes", {})
    for mode, threshold in expected_thresholds.items():
        actual = modes.get(mode, {}).get("threshold")
        if actual != threshold:
            errors.append(f"modes.{mode}.threshold must be {threshold}")

    for key in ["decision_doc", "canonical_summary", "window_manifest", "planning_script"]:
        relpath = status.get(key)
        if not relpath:
            errors.append(f"missing status.{key}")
        else:
            require_path(relpath, errors)

    bridge = profile.get("experimental_phase9", {}).get("bridge_forest", {})
    for key in ["apply_script", "training_script", "metrics", "summary"]:
        relpath = bridge.get(key)
        if not relpath:
            errors.append(f"missing experimental_phase9.bridge_forest.{key}")
        else:
            require_path(relpath, errors)

    metrics_path = ROOT / bridge.get("metrics", "")
    if metrics_path.exists():
        with metrics_path.open() as fh:
            metrics = json.load(fh)
        for key in ["rows", "positives", "negatives", "cv_roc_auc", "selected_feature_count"]:
            if key not in metrics:
                errors.append(f"metrics missing key: {key}")

    summary_path = ROOT / bridge.get("summary", "")
    if summary_path.exists():
        with summary_path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        if not rows:
            errors.append("bridge summary TSV has no data rows")
        for key in ["config", "n_windows", "total_accepted_edges", "total_hamming_errors"]:
            if rows and key not in rows[0]:
                errors.append(f"bridge summary missing column: {key}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("current Phase9D profile ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
