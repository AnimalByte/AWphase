#!/usr/bin/env python3
import json
from collections import Counter
import pysam

target_json = "data/derived/hg002_chr20_variants.real.json"
panel_bcf = "data/panels/hgdp_1kg_hg38/subset/panel_chr20_target_sites.bcf"
out_json = "results/qc/panel_conformity_audit.json"

with open(target_json) as f:
    target = json.load(f)

target_by_pos = {
    int(v["pos"]): (str(v["ref_allele"]), str(v["alt_allele"]))
    for v in target
}

bcf = pysam.VariantFile(panel_bcf)

panel_by_pos = {}
for rec in bcf.fetch():
    if rec.alts is None or len(rec.alts) != 1:
        continue
    panel_by_pos[int(rec.pos)] = (str(rec.ref), str(rec.alts[0]))

target_positions = set(target_by_pos)
panel_positions = set(panel_by_pos)

shared = sorted(target_positions & panel_positions)
target_only = sorted(target_positions - panel_positions)
panel_only = sorted(panel_positions - target_positions)

refalt_match = 0
refalt_mismatch = 0
mismatch_examples = []

for pos in shared:
    t = target_by_pos[pos]
    p = panel_by_pos[pos]
    if t == p:
        refalt_match += 1
    else:
        refalt_mismatch += 1
        if len(mismatch_examples) < 20:
            mismatch_examples.append({
                "pos": pos,
                "target_ref": t[0],
                "target_alt": t[1],
                "panel_ref": p[0],
                "panel_alt": p[1],
            })

report = {
    "n_target_sites": len(target_positions),
    "n_panel_sites": len(panel_positions),
    "n_shared_positions": len(shared),
    "n_target_only_positions": len(target_only),
    "n_panel_only_positions": len(panel_only),
    "n_refalt_matches_on_shared": refalt_match,
    "n_refalt_mismatches_on_shared": refalt_mismatch,
    "target_only_examples": target_only[:20],
    "panel_only_examples": panel_only[:20],
    "refalt_mismatch_examples": mismatch_examples,
}

with open(out_json, "w") as f:
    json.dump(report, f, indent=2)

print(json.dumps(report, indent=2))
