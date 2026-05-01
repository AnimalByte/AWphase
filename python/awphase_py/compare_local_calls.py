import json
from pathlib import Path

BASE = Path("results/runs/HG002")

def load_json(path: Path):
    with open(path) as f:
        return json.load(f)

on_local = load_json(BASE / "ancestry_on" / "local_decode_result.json")
off_local = load_json(BASE / "ancestry_off" / "local_decode_result.json")

on_calls = {c["pos"]: c for c in on_local["calls"]}
off_calls = {c["pos"]: c for c in off_local["calls"]}

all_pos = sorted(set(on_calls) | set(off_calls))

lines = []
lines.append("pos\ton_phase\toff_phase\ton_conf\toff_conf\tchanged")

for pos in all_pos:
    on = on_calls.get(pos)
    off = off_calls.get(pos)

    on_phase = on["phase_state"] if on else "NA"
    off_phase = off["phase_state"] if off else "NA"
    on_conf = f'{on["confidence"]:.6f}' if on else "NA"
    off_conf = f'{off["confidence"]:.6f}' if off else "NA"
    changed = "yes" if on_phase != off_phase else "no"

    lines.append(f"{pos}\t{on_phase}\t{off_phase}\t{on_conf}\t{off_conf}\t{changed}")

out = "\n".join(lines)
print(out)
