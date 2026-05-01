#!/usr/bin/env python3
import argparse
import csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local-calls-tsv", required=True)
    ap.add_argument("--scored-sites-tsv", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    score_by_pos = {}
    with open(args.scored_sites_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            score_by_pos[int(row["pos"])] = float(row["xgb_safe_prob"])

    rows = []
    with open(args.local_calls_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames + ["xgb_safe_prob", "xgb_gate_applied"]
        for row in reader:
            pos = int(row["pos"])
            p = score_by_pos.get(pos, 0.0)
            gate_applied = 0
            if row["phase_state"] != "0" and p < args.threshold:
                row["phase_state"] = "0"
                if "local_phase_state" in row:
                    row["local_phase_state"] = "0"
                gate_applied = 1
            row["xgb_safe_prob"] = f"{p:.6f}"
            row["xgb_gate_applied"] = gate_applied
            rows.append(row)

    with open(args.out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.out_tsv}")

if __name__ == "__main__":
    main()
