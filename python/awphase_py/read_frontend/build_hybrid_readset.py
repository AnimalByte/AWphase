#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

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

def first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"Could not find any of columns {candidates} in header: {list(cols)}")

def normalize_hp(x):
    s = str(x).strip().upper()
    if s in {"1", "HP1"}:
        return "1"
    if s in {"2", "HP2"}:
        return "2"
    return "UNK"

def as_int(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default

def as_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-readobs-json", required=True)
    ap.add_argument("--selected-reads-tsv", required=True)
    ap.add_argument("--haplotag-tsv", required=True)
    ap.add_argument("--superreads-tsv", required=True)
    ap.add_argument("--out-hybrid-json", required=True)
    ap.add_argument("--out-summary-json", required=True)
    ap.add_argument("--tagged-multiplier", type=float, default=1.20)
    ap.add_argument("--superread-base-weight", type=float, default=0.90)
    ap.add_argument("--superread-weight-step", type=float, default=0.08)
    ap.add_argument("--superread-weight-cap", type=float, default=1.35)
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

    base_obs = []
    tagged_obs_count = 0

    for o in raw_obs:
        rid = str(o.get("read_id"))
        if rid not in selected_ids:
            continue

        obs = dict(o)

        # Ensure required minimal fields exist
        obs["read_id"] = rid
        obs["site_pos"] = as_int(obs.get("site_pos"))
        obs["allele"] = as_int(obs.get("allele"))
        obs["baseq"] = as_int(obs.get("baseq"), 30)
        obs["mapq"] = as_int(obs.get("mapq"), 60)

        # Keep V2 as default backbone weight
        if "weight_v2" not in obs:
            obs["weight_v2"] = min(1.0, (obs["baseq"] / 40.0) * (obs["mapq"] / 60.0))

        # If this read is confidently tagged, modestly boost it
        if rid in tagged_ids:
            tagged_obs_count += 1
            w = as_float(obs.get("weight_v2"), 0.5) * args.tagged_multiplier
            obs["weight_v2"] = min(1.50, w)
            if "allele_confidence" not in obs:
                obs["allele_confidence"] = 0.25

        base_obs.append(obs)

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
            pos = as_int(r[sr_pos_col])
            allele = as_int(r[sr_allele_col])
        except Exception:
            continue

        support = 3
        if sr_support_col is not None:
            support = max(1, as_int(r.get(sr_support_col, 3), 3))

        weight = min(args.superread_weight_cap, args.superread_base_weight + args.superread_weight_step * support)
        conf = min(1.0, 0.20 * support)

        obs = {
            "read_id": f"HYBRID_SR::{r[sr_id_col]}",
            "site_pos": pos,
            "allele": allele,
            "baseq": 40,
            "mapq": 60,
            "raw_baseq": 40,
            "weight_v2": weight,
            "allele_confidence": conf,
            "is_ambiguous": False,
            "context_bases_compared": 0,
            "variant_class": "snv",
            "penalty_multiplier": 1.0,
        }
        super_obs.append(obs)

    hybrid_obs = base_obs + super_obs
    write_obs_json(args.out_hybrid_json, hybrid_obs, shape, key)

    summary = {
        "selected_read_ids": len(selected_ids),
        "tagged_read_ids": len(tagged_ids),
        "base_obs": len(base_obs),
        "tagged_obs_boosted": tagged_obs_count,
        "superread_obs": len(super_obs),
        "hybrid_obs_total": len(hybrid_obs),
        "tagged_multiplier": args.tagged_multiplier,
        "superread_base_weight": args.superread_base_weight,
        "superread_weight_step": args.superread_weight_step,
        "superread_weight_cap": args.superread_weight_cap,
    }
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
