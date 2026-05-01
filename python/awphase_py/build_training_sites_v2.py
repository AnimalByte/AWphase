#!/usr/bin/env python3
import argparse
import csv
import sys
import json
from pathlib import Path

def load_local_debug(path):
    with open(path) as fh:
        obj = json.load(fh)
    if isinstance(obj, dict) and "sites" in obj and isinstance(obj["sites"], list):
        return {int(s["pos"]): s for s in obj["sites"]}
    raise ValueError(f"Unsupported local debug JSON structure: {path}")

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

def load_tsv_by_pos(path):
    out = {}
    if not path:
        return out
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos = None
            for key in ("pos", "site_pos", "position", "variant_pos"):
                if key in row and str(row[key]).strip():
                    pos = int(str(row[key]).strip())
                    break
            if pos is not None:
                out[pos] = row
    return out

def get(row, keys, default=""):
    for k in keys:
        if k in row and str(row[k]).strip():
            return row[k]
    return default

def truth_label(row):
    pred = get(row, ["pred_phase_state", "predicted_phase_state", "phase_state_pred", "pred_state", "local_phase_state"])
    truth = get(row, ["truth_phase_state", "phase_state_truth", "truth_state"])
    if pred == "" or truth == "":
        return None
    try:
        pred_i = int(pred)
        truth_i = int(truth)
    except ValueError:
        return None
    if pred_i == 0:
        return 0
    return 1 if pred_i == truth_i else -1

def to_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def to_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    rows_out = []

    with open(args.manifest_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for m in reader:
            exp = m["experiment_id"]
            dbg = load_local_debug(m["local_debug_json"])
            truth = load_tsv_by_pos(m["truth_site_comparison_tsv"])
            site_dbg = load_tsv_by_pos(m.get("site_debug_tsv", ""))

            positions = sorted(set(dbg) & set(truth))
            for pos in positions:
                s = dbg[pos]
                t = truth[pos]
                sd = site_dbg.get(pos, {})

                label = truth_label(t)
                if label is None:
                    continue

                window_index = to_int(s.get("window_index", -1))
                group_id = f"{exp}__w{window_index}"

                rows_out.append({
                    "experiment_id": exp,
                    "group_id": group_id,
                    "window_index": window_index,
                    "pos": pos,
                    "label": label,
                    "phase_state": to_int(s.get("phase_state", 0)),
                    "confidence": to_float(s.get("confidence", 0.0)),
                    "alt_support": to_float(s.get("alt_support", 0.0)),
                    "ref_support": to_float(s.get("ref_support", 0.0)),
                    "read_bias": to_float(s.get("read_bias", 0.0)),
                    "donor_margin": to_float(s.get("donor_margin", 0.0)),
                    "combined_donor_bias": to_float(s.get("combined_donor_bias", 0.0)),
                    "path_sign": to_int(s.get("path_sign", 0)),
                    "read_sign": to_int(s.get("read_sign", 0)),
                    "donor_sign": to_int(s.get("donor_sign", 0)),
                    "path_margin_agrees": int(bool(s.get("path_margin_agrees", False))),
                    "donor_read_conflict": int(bool(s.get("donor_read_conflict", False))),
                    "weak_read_evidence": int(bool(s.get("weak_read_evidence", False))),
                    "forced_abstain": int(bool(s.get("forced_abstain", False))),
                    "decision_source": s.get("decision_source", ""),
                    "neighbor_dist": to_int(sd.get("neighbor_dist", 0)),
                    "reads_seen": to_int(sd.get("reads_seen", 0)),
                    "ref_obs": to_int(sd.get("ref_obs", 0)),
                    "alt_obs": to_int(sd.get("alt_obs", 0)),
                    "other_base": to_int(sd.get("other_base", 0)),
                    "deletions": to_int(sd.get("deletions", 0)),
                    "refskips": to_int(sd.get("refskips", 0)),
                    "low_mapq": to_int(sd.get("low_mapq", 0)),
                    "low_baseq": to_int(sd.get("low_baseq", 0)),
                    "secondary_or_supplementary": to_int(sd.get("secondary_or_supplementary", 0)),
                    "duplicates": to_int(sd.get("duplicates", 0)),
                    "emitted_obs": to_int(sd.get("emitted_obs", 0)),
                    "ambiguous_obs": to_int(sd.get("ambiguous_obs", 0)),
                    "alt_like_context_obs": to_int(sd.get("alt_like_context_obs", 0)),
                    "ref_like_context_obs": to_int(sd.get("ref_like_context_obs", 0)),
                    "mean_allele_confidence": to_float(sd.get("mean_allele_confidence", 0.0)),
                    "mean_signed_score_delta": to_float(sd.get("mean_signed_score_delta", 0.0)),
                    "mean_abs_score_delta": to_float(sd.get("mean_abs_score_delta", 0.0)),
                    "mean_local_ref_score": to_float(sd.get("mean_local_ref_score", 0.0)),
                    "mean_local_alt_score": to_float(sd.get("mean_local_alt_score", 0.0)),
                    "usable_site": int(str(sd.get("usable", "True")) == "True"),
                })

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows_out[0].keys())
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} rows to {out}")

if __name__ == "__main__":
    main()
