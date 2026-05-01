# Scripts

The scripts directory mixes current runners, validators, and historical
experiment helpers. Prefer these entry points for new local work:

- `phase9d_window.py`: manifest-driven Phase9D planning and local-input checks
- `validate_current_phase9d_profile.py`: current default/profile validation
- `validate_benchmark_bundles.py`: benchmark truth VCF/BED bundle validation
- `run_window_awphase_and_whatshap.sh`: heavy single-window benchmark runner
- `run_phase7a_window.sh`: Phase7A prerequisite runner used by later phases

The `phase8/` and `phase9/` subdirectories contain retained experiment scripts.
Phase9 bridge/forest scripts are experimental until a newer decision note
promotes them over the documented Phase8F threshold `0.80` default.

Download-oriented scripts are retained for provenance, but source validation and
the Phase9D planning wrapper do not invoke them.
