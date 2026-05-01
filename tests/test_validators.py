import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_benchmark_version_detects_known_bundles():
    validator = load_script("validate_benchmark_bundles")

    assert validator.benchmark_version("HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz") == "v5.0q"
    assert validator.benchmark_version("hg002_chr22_full.v5q.variants.real.json") == "v5.0q"
    assert validator.benchmark_version("HG002_GRCh38_v4.2.1_benchmark.chr22.bed") == "v4.2.1"
    assert validator.benchmark_version("truth/unknown.bed") is None


def test_manifest_pair_validator_reports_mixed_bundle(tmp_path, monkeypatch):
    validator = load_script("validate_benchmark_bundles")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "label\ttruth_vcf\ttruth_bed\n"
        "bad\tdata/HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz\t"
        "data/HG002_GRCh38_v4.2.1_benchmark.chr22.bed\n"
    )
    monkeypatch.setattr(validator, "TSV_PATHS", [manifest])

    errors = []
    validator.validate_manifest_pairs(errors)

    assert len(errors) == 1
    assert "bad mixes v5.0q VCF with v4.2.1 BED" in errors[0]


def test_current_phase9d_profile_validator_passes():
    validator = load_script("validate_current_phase9d_profile")

    assert validator.main() == 0
