#!/usr/bin/env python3
import argparse
import json

def load_positions(path):
    keep = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            keep.add(int(parts[1]))
    return keep

def filter_obj(obj, keep):
    # list of site dicts
    if isinstance(obj, list):
        if len(obj) == 0:
            return obj
        if isinstance(obj[0], dict) and "pos" in obj[0]:
            return [x for x in obj if int(x["pos"]) in keep]
        return obj

    # dict wrapper with site list
    if isinstance(obj, dict):
        for key in ("sites", "priors", "records", "items"):
            if key in obj and isinstance(obj[key], list):
                first = obj[key][0] if obj[key] else None
                if isinstance(first, dict) and "pos" in first:
                    out = dict(obj)
                    out[key] = [x for x in obj[key] if int(x["pos"]) in keep]
                    return out

        # dict keyed by positions
        if all(str(k).isdigit() for k in obj.keys()):
            return {k: v for k, v in obj.items() if int(k) in keep}

        # dict keyed by chrom -> nested structure
        out = {}
        changed = False
        for k, v in obj.items():
            fv = filter_obj(v, keep)
            out[k] = fv
            if fv is not v:
                changed = True
        if changed:
            return out

    raise ValueError(
        f"Unsupported site priors JSON structure: top-level type={type(obj).__name__}"
    )

def count_sites(obj):
    if isinstance(obj, list):
        if len(obj) == 0:
            return 0
        if isinstance(obj[0], dict) and "pos" in obj[0]:
            return len(obj)
        return 0
    if isinstance(obj, dict):
        for key in ("sites", "priors", "records", "items"):
            if key in obj and isinstance(obj[key], list):
                first = obj[key][0] if obj[key] else None
                if first is None or (isinstance(first, dict) and "pos" in first):
                    return len(obj[key])
        if all(str(k).isdigit() for k in obj.keys()):
            return len(obj)
        total = 0
        for v in obj.values():
            total += count_sites(v)
        return total
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-priors-json", required=True)
    ap.add_argument("--positions-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    keep = load_positions(args.positions_tsv)

    with open(args.site_priors_json) as fh:
        obj = json.load(fh)

    out = filter_obj(obj, keep)

    with open(args.out_json, "w") as fh:
        json.dump(out, fh)

    print(f"filtered_site_priors={count_sites(out)}")

if __name__ == "__main__":
    main()
