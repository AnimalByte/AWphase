#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

def walk(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for item in obj:
            yield from walk(item)

def find_candidate_container(obj):
    if isinstance(obj, dict) and ("donors" in obj or "composites" in obj):
        return obj
    for d in walk(obj):
        if isinstance(d, dict) and ("donors" in d or "composites" in d):
            return d
    return {"donors": [], "composites": []}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    obj = load_json(args.json)
    container = find_candidate_container(obj)
    donors = container.get("donors", []) or []
    composites = container.get("composites", []) or []

    rows = []

    for rank, d in enumerate(donors, start=1):
        rows.append({
            "kind": "donor",
            "rank": rank,
            "donor_id": d.get("donor_id", ""),
            "hap_id": d.get("hap_id", ""),
            "group": d.get("group", ""),
            "raw_score": d.get("raw_score", ""),
            "ancestry_weighted_score": d.get("ancestry_weighted_score", ""),
            "matched_sites": d.get("matched_sites", ""),
            "total_sites": d.get("total_sites", ""),
            "state_id": "",
            "members": "",
            "score": "",
        })

    for rank, c in enumerate(composites, start=1):
        members = c.get("members", [])
        member_str = ",".join(
            f"{m[0]}:{m[1]}" if isinstance(m, list) and len(m) == 2 else str(m)
            for m in members
        ) if isinstance(members, list) else ""
        rows.append({
            "kind": "composite",
            "rank": rank,
            "donor_id": "",
            "hap_id": "",
            "group": "",
            "raw_score": "",
            "ancestry_weighted_score": "",
            "matched_sites": "",
            "total_sites": "",
            "state_id": c.get("state_id", ""),
            "members": member_str,
            "score": c.get("score", ""),
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "kind",
                "rank",
                "donor_id",
                "hap_id",
                "group",
                "raw_score",
                "ancestry_weighted_score",
                "matched_sites",
                "total_sites",
                "state_id",
                "members",
                "score",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out}")

if __name__ == "__main__":
    main()
