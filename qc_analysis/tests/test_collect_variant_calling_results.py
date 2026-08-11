import csv
import gzip
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "qc_analysis/scripts/collect_variant_calling_results.py"
WRAPPER = ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"


def prepare(tmp_path, missing=(), single=False):
    source = tmp_path / "aggregate"
    for name in ("vcf", "mtcn", "round2_coverage", "numt_decoy_coverage"):
        (source / name).mkdir(parents=True)
    metadata = tmp_path / "samples.tsv"
    metadata.write_text("sample\tspecies\nS1\tSpecies one\n")
    out = tmp_path / "out"
    config = tmp_path / "qc.yaml"
    config.write_text(
        "collect_variant_calling:\n"
        f"  input_root: {source}\n  outdir: {out}\n"
        "  vcf_subdir: vcf\n  mtcn_subdir: mtcn\n"
        "  round2_coverage_subdir: round2_coverage\n"
        "  numt_decoy_coverage_subdir: numt_decoy_coverage\n"
        f"  metadata: {metadata}\n  allow_single_cov: {'true' if single else 'false'}\n"
    )
    (source / "mtcn/S1.round2.mtcn.tsv").write_text(
        "mt_median_coverage\tnuclear_median_coverage\tmtcn_median\n300\t20\t30\n"
    )
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchrM\t1\t.\tA\tG\t.\tPASS\t.\tDP:AF\t101:0.1\n"
    (source / "vcf/S1.round2.original_coords.clean.final.split.vcf").write_text(vcf)
    with gzip.open(source / "vcf/S1.round2.original_coords.clean.final.split.vcf.gz", "wt") as handle:
        handle.write(vcf)
    if "round2" not in missing:
        (source / "round2_coverage/S1.round2.original_coords.per_base_coverage.tsv").write_text(
            "chrM\t1\tchrM\t100\nchrM\t2\tchrM\t500\nchrM\t3\tchrM\t300\nchrM\t4\tchrM\t50\n"
        )
    if "numt" not in missing:
        (source / "numt_decoy_coverage/S1.numt_decoy.clean.realigned.per_base_coverage.tsv").write_text(
            "chrM\t1\tchrM\t200\nchrM\t2\tchrM\t400\nchrM\t3\tchrM\t300\nchrM\t5\tchrM\t75\n"
        )
    return config, out


def run(config, wrapper=False):
    cmd = (["bash", str(WRAPPER), "collect_variant_calling_results", str(config)] if wrapper
           else ["python3", str(SCRIPT), "--config", str(config)])
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)


def summary(out):
    with (out / "reports/variant_calling_collection_summary.tsv").open() as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_aggregated_merge_retains_unique_coordinates_and_prefers_gzip(tmp_path):
    config, out = prepare(tmp_path)
    result = run(config)
    assert result.returncode == 0, result.stderr
    assert summary(out)["status"] == "OK"
    assert (out / "collected_vcf/S1.round2.original_coords.clean.final.split.vcf.gz").is_symlink()
    rows = list(csv.reader((out / "collected_cov/S1.merged.max_depth.per_base_coverage.tsv").open(), delimiter="\t"))
    assert rows[1:] == [["chrM", "1", "chrM", "200"], ["chrM", "2", "chrM", "500"],
                       ["chrM", "3", "chrM", "300"], ["chrM", "4", "chrM", "50"],
                       ["chrM", "5", "chrM", "75"]]
    assert "merged_coverage_written=1" in result.stdout


def test_each_missing_coverage_source_is_unambiguous(tmp_path):
    for missing, label in (("round2", "round2_coverage"), ("numt", "numt_decoy_coverage")):
        config, out = prepare(tmp_path / missing, missing=(missing,))
        assert run(config).returncode == 0
        row = summary(out)
        assert row["status"] == "MISSING_INPUT" and row["missing_files"] == label
        assert not (out / "collected_cov/S1.merged.max_depth.per_base_coverage.tsv").exists()


def test_allow_single_coverage_writes_output(tmp_path):
    config, out = prepare(tmp_path, missing=("round2",), single=True)
    assert run(config).returncode == 0
    row = summary(out)
    assert row["status"] == "OK" and row["missing_files"] == "NA"
    assert (out / "collected_cov/S1.merged.max_depth.per_base_coverage.tsv").exists()


def test_wrapper_direct_mode_and_standardized_liftover_names(tmp_path):
    config, out = prepare(tmp_path)
    result = run(config, wrapper=True)
    assert result.returncode == 0, result.stderr
    assert (out / "collected_vcf/S1.round2.original_coords.clean.final.split.vcf.gz").exists()
    assert (out / "collected_cov/S1.merged.max_depth.per_base_coverage.tsv").exists()


def test_all_aggregate_directories_are_validated(tmp_path):
    config, _ = prepare(tmp_path)
    (tmp_path / "aggregate/mtcn/S1.round2.mtcn.tsv").unlink()
    (tmp_path / "aggregate/mtcn").rmdir()
    result = run(config)
    assert result.returncode != 0
    assert "mtcn=" in result.stderr
