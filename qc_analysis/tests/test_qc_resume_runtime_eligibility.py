import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "qc_analysis/scripts/qc_array_manifest.py"
WRAPPER = ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"
ELIGIBILITY = ROOT / "qc_analysis/scripts/qc_sample_runtime_eligibility.py"
SAMPLE_QC = ROOT / "qc_analysis/scripts/run_sample_variant_filtering.py"


def parse_values(stdout: str) -> dict[str, str]:
    values = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def test_completed_liftover_is_revalidated_without_mt_codon_tag(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tspecies\nS1\tSpecies_one\n")
    out = tmp_path / "liftover"
    lifted = out / "vcf_lifted_raw/S1.lifted.raw.vcf"
    lifted.parent.mkdir(parents=True)
    lifted.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=SRC_POS,Number=1,Type=Integer,Description=\"Original species position\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    config = tmp_path / "qc.yaml"
    config.write_text(
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    sample_ref_file: {samples}\n"
        f"    output_dir: {out}\n"
    )

    result = subprocess.run(
        [sys.executable, str(MANIFEST), "coordinate_liftover", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    values = parse_values(result.stdout)
    assert values["COUNT"] == "1"
    assert values["STATE"] == "scheduled"
    assert Path(values["TASK_FILE"]).read_text() == "S1\n"
    assert "runtime_revalidate" in Path(values["MANIFEST"]).read_text()


def test_missing_sample_candidates_are_not_silently_treated_as_complete(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tspecies\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    sample_ref_file: {samples}\n"
        f"    output_dir: {tmp_path / 'liftover'}\n"
    )
    result = subprocess.run(
        [sys.executable, str(MANIFEST), "coordinate_liftover", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode != 0
    assert "no candidate samples for coordinate_liftover" in result.stderr


def test_sample_qc_recomputes_and_replaces_existing_report(tmp_path):
    collection = tmp_path / "collection.tsv"
    collection.write_text(
        "sample\tspecies\tmt_median_coverage\tPercent_100\tnuclear_median_coverage\tmtcn_median\tMAD\tstatus\n"
        "S1\tSpecies_one\t150\t95\t30\t60\t0.1\tOK\n"
    )
    out = tmp_path / "sample_qc"
    report = out / "reports/sample_qc.tsv"
    report.parent.mkdir(parents=True)
    report.write_text("sample\tqc_status\nSTALE\tFAIL\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "sample_variant_filtering:\n"
        "  enabled: true\n"
        f"  input_summary: {collection}\n"
        f"  output_dir: {out}\n"
        "  thresholds:\n"
        "    mt_median_coverage_min: 100\n"
        "    percent_100_min: 90\n"
        "    nuclear_median_coverage_min: 20\n"
        "    mtcn_min: 40\n"
        "    mad_max: 0.5\n"
    )
    result = subprocess.run(
        [sys.executable, str(SAMPLE_QC), "--config", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert "replacing existing report" in result.stderr
    with report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample"] for row in rows] == ["S1"]
    assert rows[0]["qc_status"] == "PASS"


def write_runtime_config(tmp_path: Path, mapping_rows: str) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    mapping = tmp_path / "codon_map.tsv"
    mapping.write_text("sample\treference_key\n" + mapping_rows)
    liftover = tmp_path / "liftover"
    inputs = liftover / "vcf_lifted_raw"
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "S1.lifted.raw.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    reports = liftover / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "S1.coordinate_liftover_qc.tsv").write_text("sample\tS1\nstatus\tcompleted\n")
    collection = tmp_path / "collection"
    (collection / "reports").mkdir(parents=True, exist_ok=True)
    (collection / "reports/variant_calling_collection_summary.tsv").write_text(
        "sample\tstatus\nS1\tOK\n"
    )
    config = tmp_path / "runtime.yaml"
    config.write_text(
        "collect_variant_calling:\n"
        f"  outdir: {collection}\n"
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    output_dir: {liftover}\n"
        "codon_match:\n"
        "  paths:\n"
        f"    sample_reference_map: {mapping}\n"
        f"    input_vcf_dir: {inputs}\n"
        f"    output_dir: {tmp_path / 'codon'}\n"
        f"    reports_dir: {tmp_path / 'codon/reports'}\n"
        "  settings:\n"
        "    input_vcf_pattern: \"{sample}.lifted.raw.vcf\"\n"
        "    output_suffix: \".lifted.codon.vcf\"\n"
    )
    return config


def test_codon_runtime_eligibility_uses_current_pass_production_map(tmp_path):
    eligible_config = write_runtime_config(tmp_path, "S1\tmtref_abc\n")
    eligible = subprocess.run(
        [sys.executable, str(ELIGIBILITY), "codon_match", "S1", str(eligible_config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert eligible.returncode == 0, eligible.stderr
    assert "ELIGIBLE=1" in eligible.stdout

    ineligible_config = write_runtime_config(tmp_path / "excluded", "OTHER\tmtref_def\n")
    ineligible = subprocess.run(
        [sys.executable, str(ELIGIBILITY), "codon_match", "S1", str(ineligible_config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert ineligible.returncode == 0, ineligible.stderr
    assert "ELIGIBLE=0" in ineligible.stdout
    assert "sample_not_in_pass_production_codon_map" in ineligible.stdout


def test_array_worker_runtime_exclusion_is_successful_skip(tmp_path):
    config = write_runtime_config(tmp_path, "OTHER\tmtref_def\n")
    tasks = tmp_path / "tasks.txt"
    tasks.write_text("S1\n")
    result = subprocess.run(
        [
            "bash", str(WRAPPER), "--array-task", "--runtime-eligibility",
            "--task-file", str(tasks), "codon_match", str(config),
        ],
        cwd=ROOT, text=True, capture_output=True,
        env={**os.environ, "SLURM_ARRAY_TASK_ID": "1"},
    )
    assert result.returncode == 0, result.stderr
    assert "runtime_status=skipped" in result.stderr
    assert "sample_not_in_pass_production_codon_map" in result.stderr


def test_runtime_exclusion_removes_stale_codon_outputs(tmp_path):
    config = write_runtime_config(tmp_path, "OTHER\tmtref_def\n")
    stale_vcf = tmp_path / "codon/vcf_codon/S1.lifted.codon.vcf"
    stale_summary = tmp_path / "codon/reports/S1.codon_match_summary.tsv"
    stale_vcf.parent.mkdir(parents=True, exist_ok=True)
    stale_summary.parent.mkdir(parents=True, exist_ok=True)
    stale_vcf.write_text("stale\n")
    stale_summary.write_text("sample\tstatus\nS1\tcompleted\n")
    result = subprocess.run(
        [sys.executable, str(ELIGIBILITY), "codon_match", "S1", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert "ELIGIBLE=0" in result.stdout
    assert not stale_vcf.exists()
    assert not stale_summary.exists()
