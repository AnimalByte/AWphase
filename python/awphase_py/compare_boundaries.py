import json
from pathlib import Path

BASE = Path("results/runs/HG002")

def load_json(path: Path):
    with open(path) as f:
        return json.load(f)

on = load_json(BASE / "ancestry_on" / "boundary_decisions.json")
off = load_json(BASE / "ancestry_off" / "boundary_decisions.json")

lines = []
lines.append("left_block\tright_block\ton_decision\toff_decision\ton_conf\toff_conf\tchanged")

for on_b, off_b in zip(on, off):
    left = on_b["left_block"]
    right = on_b["right_block"]
    on_dec = on_b["decision"]
    off_dec = off_b["decision"]
    on_conf = f'{on_b["confidence"]["value"]:.6f}'
    off_conf = f'{off_b["confidence"]["value"]:.6f}'
    changed = "yes" if on_dec != off_dec else "no"

    lines.append(f"{left}\t{right}\t{on_dec}\t{off_dec}\t{on_conf}\t{off_conf}\t{changed}")

out = "\n".join(lines)
print(out)
