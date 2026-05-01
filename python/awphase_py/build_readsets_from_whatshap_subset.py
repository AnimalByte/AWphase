#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"Could not find any of columns {candidates} in header: {list(cols)}")

def load_obs_json(path):
    obj = json.load(open(path))
    if isinstance(obj, list):
        return obj, "list", None
    if isinstance(obj, dict):
        for key in ["observations", "read_observations", "reads"]:
            if key in obj and isinstance(obj[key], list):
                return obj[key], "dict", key
    raise SystemExit(f"Unsupported readobs JSON shape in {path}")

def write_obs_json(path, obs, shape, key):
    if shape == "list":
        json.dump(obs, open(path, "w"), indent=2)
    else:
        json.dump({key: obs}, open(path, "w"), indent=2)

def normalize_hp(x):
    s = str(x).strip().upper()
    if s in {"1", "HP1"}:
        return "1"
    if s in {"2", "HP2"}:
        return "2"
    return "UNK"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-readobs-json", required=True)
    ap.add_argument("--selected-reads-tsv", required=True)
    ap.add_argument("--haplotag-tsv", required=True)
    ap.add_argument("--superreads-tsv", required=True)
    ap.add_argument("--out-selected-json", required=True)
    ap.add_argument("--out-tagged-json", required=True)
    ap.add_argument("--out-superreads-json", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    raw_obs, shape, key = load_obs_json(args.raw_readobs_json)

    selected_rows = read_tsv(args.selected_reads_tsv)
    haplo_rows = read_tsv(args.haplotag_tsv)
    super_rows = read_tsv(args.superreads_tsv)

    sel_read_col = first_present(selected_rows[0].keys(), ["read_id", "read", "name"])
    hp_read_col  = first_present(haplo_rows[0].keys(), ["read_id", "read", "name"])
    hp_col       = first_present(haplo_rows[0].keys(), ["hp", "HP", "haplotype"])

    selected_ids = {r[sel_read_col] for r in selected_rows}
    tagged_ids = {
        r[hp_read_col]
        for r in haplo_rows
        if normalize_hp(r[hp_col]) in {"1", "2"}
    }

    selected_obs = [o for o in raw_obs if str(o.get("read_id")) in selected_ids]
    tagged_obs   = [o for o in raw_obs if str(o.get("read_id")) in tagged_ids]

    sr_cols = super_rows[0].keys()
    sr_id_col = first_present(sr_cols, ["superread_id", "read_id", "sr_id", "name"])
    sr_pos_col = first_present(sr_cols, ["site_pos", "pos"])
    sr_allele_col = first_present(sr_cols, ["allele"])
    try:
        sr_support_col = first_present(sr_cols, ["support_reads", "n_reads", "support", "read_count", "count"])
    except SystemExit:
        sr_support_col = None

    super_obs = []
    for r in super_rows:
        try:
            pos = int(r[sr_pos_col])
            allele = int(r[sr_allele_col])
        except Exception:
            continue
        support = 3
        if sr_support_col is not None:
            try:
                support = max(1, int(float(r[sr_support_col])))
            except Exception:
                support = 3

        obs = {
            "read_id": f"SR::{r[sr_id_col]}",
            "site_pos": pos,
            "allele": allele,
            "baseq": 40,
            "mapq": 60,
            "raw_baseq": 40,
            "weight_v2": min(1.25, 0.50 + 0.10 * support),
            "allele_confidence": min(1.0, 0.20 * support),
            "is_ambiguous": False,
            "context_bases_compared": 0,
            "variant_class": "snv",
            "penalty_multiplier": 1.0,
        }
        super_obs.append(obs)

    write_obs_json(args.out_selected_json, selected_obs, shape, key)
    write_obs_json(args.out_tagged_json, tagged_obs, shape, key)
    write_obs_json(args.out_superreads_json, super_obs, shape, key)

    summary = {
        "raw_obs": len(raw_obs),
        "selected_obs": len(selected_obs),
        "tagged_obs": len(tagged_obs),
        "superread_obs": len(super_obs),
        "selected_read_ids": len(selected_ids),
        "tagged_read_ids": len(tagged_ids),
        "distinct_superreads": len({r[sr_id_col] for r in super_rows}),
    }
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
