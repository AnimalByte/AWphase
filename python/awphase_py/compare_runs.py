import json
from pathlib import Path

BASE = Path("results/runs/HG002")

def load_json(path: Path):
    with open(path) as f:
        return json.load(f)

on_summary = load_json(BASE / "ancestry_on" / "benchmark_summary.json")
off_summary = load_json(BASE / "ancestry_off" / "benchmark_summary.json")

on_boundary = load_json(BASE / "ancestry_on" / "boundary_decision.json")
off_boundary = load_json(BASE / "ancestry_off" / "boundary_decision.json")

print("awphase run comparison")
print("======================")
print(f"sample: {on_summary['sample_id']}")
print(f"region: {on_summary['region']}")
print()

print("local summary")
print("-------------")
print(f"ancestry_on : phased={on_summary['number_phased']} abstained={on_summary['number_abstained']}")
print(f"ancestry_off: phased={off_summary['number_phased']} abstained={off_summary['number_abstained']}")
print()

print("bridge summary")
print("--------------")
print(f"ancestry_on : decision={on_boundary['decision']} confidence={on_boundary['confidence']['value']:.6f}")
print(f"ancestry_off: decision={off_boundary['decision']} confidence={off_boundary['confidence']['value']:.6f}")
print()

delta_phased = on_summary["number_phased"] - off_summary["number_phased"]
delta_abstained = on_summary["number_abstained"] - off_summary["number_abstained"]

print("difference")
print("----------")
print(f"delta_phased: {delta_phased:+d}")
print(f"delta_abstained: {delta_abstained:+d}")
