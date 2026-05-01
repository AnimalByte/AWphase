#!/usr/bin/env python3
import argparse, csv, json
from collections import defaultdict
from pathlib import Path

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def as_int(x, d=0):
    try: return int(float(x))
    except Exception: return d

def as_float(x, d=0.0):
    try: return float(x)
    except Exception: return d

def detect_score_col(row):
    for c in ("rank_score","pred_score","score","candidate_score","xgb_score","prob"):
        if c in row: return c
    for c in row:
        if "score" in c.lower():
            return c
    raise SystemExit(f"Could not find score column in candidate scores: {list(row.keys())}")

def load_regions(path):
    out = {}
    for r in read_tsv(path):
        out[r["region_id"]] = {
            "positions": [as_int(x) for x in r["positions"].split(",") if str(x).strip()],
            "left_anchor_state": as_int(r.get("left_anchor_state",0)),
            "right_anchor_state": as_int(r.get("right_anchor_state",0)),
        }
    return out

def load_candidate_best(path):
    rows = read_tsv(path)
    score_col = detect_score_col(rows[0])
    by_pos = defaultdict(list)
    for r in rows:
        pos = as_int(r.get("pos", r.get("site_pos")))
        by_pos[pos].append({
            "state": as_int(r.get("candidate_phase_state",0)),
            "score": as_float(r.get(score_col,0.0)),
            "donor_bias": as_float(r.get("f_donor_bias", r.get("donor_bias",0.0))),
            "anchor_agreement": as_float(r.get("f_anchor_agreement", r.get("anchor_agreement",0.0))),
            "local_state": as_int(r.get("f_candidate_matches_local",0)),  # binary proxy
        })
    out = {}
    for pos,lst in by_pos.items():
        lst = sorted(lst, key=lambda x:x["score"], reverse=True)
        second = lst[1]["score"] if len(lst)>1 else -999.0
        out[pos] = {
            "top_state": lst[0]["state"],
            "margin": lst[0]["score"] - second,
            "donor_bias": lst[0]["donor_bias"],
            "anchor_agreement": lst[0]["anchor_agreement"],
        }
    return out

def load_fragments(path):
    by_region = defaultdict(list)
    for r in read_tsv(path):
        by_region[r["region_id"]].append(r)
    return by_region

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions-tsv", required=True)
    ap.add_argument("--active-reads-tsv", required=True)
    ap.add_argument("--candidate-scores-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary-json", required=True)
    args = ap.parse_args()

    regions = load_regions(args.regions_tsv)
    cands = load_candidate_best(args.candidate_scores_tsv)
    frags = load_fragments(args.active_reads_tsv)

    rows = []
    for rid, reg in regions.items():
        pos_meta = []
        for p in reg["positions"]:
            meta = cands.get(p, {"top_state":0,"margin":0.0,"donor_bias":0.0,"anchor_agreement":0.0})
            pos_meta.append({
                "pos": p,
                "top_state": meta["top_state"],
                "margin": meta["margin"],
                "donor_bias": meta["donor_bias"],
                "anchor_agreement": meta["anchor_agreement"],
            })
        rows.append({
            "region_id": rid,
            "n_sites": len(reg["positions"]),
            "positions_json": json.dumps(reg["positions"]),
            "left_anchor_state": reg["left_anchor_state"],
            "right_anchor_state": reg["right_anchor_state"],
            "site_meta_json": json.dumps(pos_meta),
            "fragments_json": json.dumps(frags.get(rid, [])),
        })

    write_tsv(args.out_tsv, rows, list(rows[0].keys()) if rows else ["region_id","positions_json"])
    summary = {
        "regions": len(rows),
        "rows": len(rows),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
