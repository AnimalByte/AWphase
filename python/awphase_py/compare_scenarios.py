import json
from pathlib import Path

runs = [
    ("clean_on", Path("results/runs/HG002_clean/ancestry_on/run_report.json")),
    ("clean_off", Path("results/runs/HG002_clean/ancestry_off/run_report.json")),
    ("noisy_on", Path("results/runs/HG002_noisy/ancestry_on/run_report.json")),
    ("noisy_off", Path("results/runs/HG002_noisy/ancestry_off/run_report.json")),
    ("messy_on", Path("results/runs/HG002_messy/ancestry_on/run_report.json")),
    ("messy_off", Path("results/runs/HG002_messy/ancestry_off/run_report.json")),
]

data = {}
for label, path in runs:
    payload = json.loads(path.read_text())
    data[label] = payload

print("awphase scenario comparison")
print("===========================")
print()

print("label\tjoined\tflipped\tabstained\tphased\tlocal_abstained")
for label in ["clean_on", "clean_off", "noisy_on", "noisy_off", "messy_on", "messy_off"]:
    s = data[label]["summary"]
    print(
        f"{label}\t{s['joined_boundaries']}\t{s['flipped_boundaries']}\t"
        f"{s['abstained_boundaries']}\t{s['number_phased']}\t{s['number_abstained']}"
    )

print()
print("block coherence")
print("---------------")
print("label\tblock_1\tblock_2\tblock_3\tblock_4")
for label in ["clean_on", "clean_off", "noisy_on", "noisy_off", "messy_on", "messy_off"]:
    blocks = data[label]["derived_block_summaries"]
    vals = [f"{b['coherence']:.3f}" for b in blocks]
    print(f"{label}\t" + "\t".join(vals))

print()
print("block informative sites")
print("-----------------------")
print("label\tblock_1\tblock_2\tblock_3\tblock_4")
for label in ["clean_on", "clean_off", "noisy_on", "noisy_off", "messy_on", "messy_off"]:
    blocks = data[label]["derived_block_summaries"]
    vals = [str(b["informative_sites"]) for b in blocks]
    print(f"{label}\t" + "\t".join(vals))
