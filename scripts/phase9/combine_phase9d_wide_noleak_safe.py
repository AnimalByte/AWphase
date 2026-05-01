#!/usr/bin/env python3
import csv
import json
from pathlib import Path

MANIFEST = Path("results/phase8f_manifests/phase8f_windows.window_beds.tsv")
OUT = Path("results/phase9/training.candidates.phase9d.block_bridges.wide.noleak.tsv")
REPORT = Path("results/phase9/training.candidates.phase9d.block_bridges.wide.noleak.report.json")

LEAK_COLS = {
    "true_relation",
    "block_a_orientation_sites",
    "block_b_orientation_sites",
    "block_a_orientation_margin",
    "block_b_orientation_margin",
}

def clean(x):
    return x.replace("\r", "") if isinstance(x, str) else x

def main():
    if not MANIFEST.exists():
        raise SystemExit(f"MISSING manifest: {MANIFEST}")

    rows = []
    fields = []
    seen = set()
    used = []
    missing = []
    empty = []
    per_file = []

    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k: clean(v) for k, v in r.items()}
            if r.get("split") != "train":
                continue

            label = r["label"]
            f = Path(f"results/phase7a_windows/{label}_illumina30x_phase7a/phase9d_block_bridge_candidates.wide.tsv")

            if not f.exists():
                missing.append(str(f))
                continue

            if f.stat().st_size == 0:
                empty.append(str(f))
                continue

            n = 0
            with open(f) as tfh:
                reader = csv.DictReader(tfh, delimiter="\t")

                if not reader.fieldnames:
                    empty.append(str(f))
                    continue

                for name in reader.fieldnames:
                    if name in LEAK_COLS:
                        continue
                    if name not in seen:
                        fields.append(name)
                        seen.add(name)

                for rr in reader:
                    rows.append({k: v for k, v in rr.items() if k not in LEAK_COLS})
                    n += 1

            per_file.append({"label": label, "path": str(f), "rows": n})
            if n > 0:
                used.append(label)

    if not rows:
        print("ERROR: no Phase9D wide candidate rows found.")
        print("You need these files first:")
        print("  results/phase7a_windows/<label>_illumina30x_phase7a/phase9d_block_bridge_candidates.wide.tsv")
        print()
        print("Missing examples:")
        for x in missing[:20]:
            print("  " + x)
        print()
        print("Empty examples:")
        for x in empty[:20]:
            print("  " + x)
        raise SystemExit(2)

    if "label" not in fields:
        raise SystemExit("ERROR: combined table would not contain a label column.")
    if "source_window" not in fields:
        raise SystemExit("ERROR: combined table would not contain source_window column.")

    OUT.parent.mkdir(parents=True, exist_ok=True)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    pos = 0
    neg = 0
    bad_label = 0

    for r in rows:
        try:
            y = int(float(r.get("label", "")))
            if y == 1:
                pos += 1
            elif y == 0:
                neg += 1
            else:
                bad_label += 1
        except Exception:
            bad_label += 1

    report = {
        "out": str(OUT),
        "rows": len(rows),
        "positives": pos,
        "negatives": neg,
        "bad_label_rows": bad_label,
        "positive_rate": pos / len(rows) if rows else None,
        "field_count": len(fields),
        "used_windows": used,
        "used_window_count": len(used),
        "missing_count": len(missing),
        "empty_count": len(empty),
        "per_file": per_file,
        "removed_leak_cols": sorted(LEAK_COLS),
    }

    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    print(json.dumps({
        "out": str(OUT),
        "rows": len(rows),
        "positives": pos,
        "negatives": neg,
        "positive_rate": report["positive_rate"],
        "used_window_count": len(used),
        "missing_count": len(missing),
        "empty_count": len(empty),
        "report": str(REPORT),
    }, indent=2))

if __name__ == "__main__":
    main()
