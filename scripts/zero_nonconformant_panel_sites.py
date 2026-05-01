#!/usr/bin/env python3
import json
from pathlib import Path
import pysam

target_json = Path("data/derived/hg002_chr20_variants.real.json")
panel_bcf = Path("data/panels/hgdp_1kg_hg38/subset/panel_chr20_target_sites.bcf")
panel_json_in = Path("data/derived/panel_chr20_candidates.json")
panel_json_out = Path("data/derived/panel_chr20_candidates.exact_conform.json")
report_out = Path("results/qc/panel_exact_conform_report.json")

with target_json.open() as f:
    target = json.load(f)

target_by_pos = {
    int(v["pos"]): (str(v["ref_allele"]), str(v["alt_allele"]))
    for v in target
}
target_positions = [int(v["pos"]) for v in target]

bcf = pysam.VariantFile(str(panel_bcf))
panel_by_pos = {}
for rec in bcf.fetch():
    if rec.alts is None or len(rec.alts) != 1:
        continue
    panel_by_pos[int(rec.pos)] = (str(rec.ref), str(rec.alts[0]))

good_idx = []
bad_idx = []
target_only_idx = []

for i, pos in enumerate(target_positions):
    t = target_by_pos[pos]
    p = panel_by_pos.get(pos)
    if p is None:
        target_only_idx.append(i)
    elif p == t:
        good_idx.append(i)
    else:
        bad_idx.append(i)

with panel_json_in.open() as f:
    panel = json.load(f)

nonzero_before = 0
nonzero_after = 0

bad_set = set(bad_idx)

for hap in panel:
    alleles = hap["alleles"]
    for a in alleles:
        if a != 0:
            nonzero_before += 1
    for i in bad_set:
        alleles[i] = 0
    for a in alleles:
        if a != 0:
            nonzero_after += 1

with panel_json_out.open("w") as f:
    json.dump(panel, f)

report = {
    "n_target_sites": len(target_positions),
    "n_exact_conformant_sites": len(good_idx),
    "n_target_only_sites": len(target_only_idx),
    "n_refalt_mismatch_sites_zeroed": len(bad_idx),
    "nonzero_panel_entries_before": nonzero_before,
    "nonzero_panel_entries_after": nonzero_after,
    "target_only_example_positions": [target_positions[i] for i in target_only_idx[:20]],
    "mismatch_example_positions": [target_positions[i] for i in bad_idx[:20]],
}
with report_out.open("w") as f:
    json.dump(report, f, indent=2)

print(json.dumps(report, indent=2))
print(f"wrote {panel_json_out}")
print(f"wrote {report_out}")
