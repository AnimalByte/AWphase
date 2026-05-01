import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def write_manifest(path: Path, truth_bed: str | None = None) -> None:
    bed = truth_bed or "data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22_20_25mb.benchmark.bed"
    path.write_text(
        "\t".join(
            [
                "split",
                "chrom",
                "start",
                "end",
                "label",
                "bam",
                "panel",
                "json",
                "truth_vcf",
                "truth_bed",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "train",
                "chr22",
                "20000000",
                "24999999",
                "chr22_20_25mb",
                "data/raw/hg002_chr22/illumina.bam",
                "data/panels/hgdp1kgp_chr22.bcf",
                "data/derived/hg002_chr22_full.v5q.variants.real.json",
                "data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz",
                bed,
            ]
        )
        + "\n"
    )


def test_load_manifest_select_and_plan(tmp_path):
    phase9d = load_script("phase9d_window")
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest)

    rows = phase9d.load_manifest(manifest)
    row = phase9d.select_window(rows, "chr22_20_25mb")
    plan = phase9d.WindowPlan(
        row=row,
        mode="default",
        threshold="0.80",
        ref_fasta="references/test.fa",
        output_root="results/phase9d_windows",
    )

    assert phase9d.sorted_labels(rows) == ["chr22_20_25mb"]
    assert plan.source_window == "chr22_20_25mb_illumina30x_phase7a"
    assert plan.phase8f_calls.endswith("/t0.80/local_calls.phase8f_xgb.tsv")
    assert "Phase8F_ancestry_weighted_multiwindow" in phase9d.render_text_plan(plan)


def test_truth_bundle_mismatch_is_reported(tmp_path):
    phase9d = load_script("phase9d_window")
    manifest = tmp_path / "manifest.tsv"
    write_manifest(
        manifest,
        truth_bed="data/truth/hg002_chr22/HG002_GRCh38_v4.2.1_benchmark.chr22.bed",
    )

    row = phase9d.load_manifest(manifest)[0]

    assert phase9d.truth_bundle_errors(row) == [
        "truth bundle version mismatch: "
        "data/truth/hg002_chr22/HG002_GRCh38_v5.0q_smvar.chr22.vcf.gz is v5.0q, "
        "data/truth/hg002_chr22/HG002_GRCh38_v4.2.1_benchmark.chr22.bed is v4.2.1"
    ]


def test_build_plan_uses_mode_threshold(tmp_path):
    phase9d = load_script("phase9d_window")
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest)

    args = SimpleNamespace(
        manifest=str(manifest),
        label="chr22_20_25mb",
        mode="conservative",
        threshold=None,
        ref_fasta="references/test.fa",
        output_root="results/phase9d_windows",
    )
    plan = phase9d.build_plan(args)

    assert plan.threshold == "0.90"


def test_check_reports_missing_external_inputs(tmp_path, capsys):
    phase9d = load_script("phase9d_window")
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest)
    args = SimpleNamespace(
        manifest=str(manifest),
        label="chr22_20_25mb",
        mode="default",
        threshold=None,
        ref_fasta="references/test.fa",
        output_root="results/phase9d_windows",
        include_phase9d_inputs=False,
    )

    assert phase9d.cmd_check(args) == 1
    captured = capsys.readouterr()
    assert "MISSING illumina_bam" in captured.err
