import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "qc_analysis/scripts/qc_array_manifest.py"
WRAPPER = ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh"
ELIGIBILITY = ROOT / "qc_analysis/scripts/qc_sample_runtime_eligibility.py"


def parse_values(stdout: str) -> dict[str, str]:
    values = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def test_completed_liftover_is_detected_without_mt_codon_tag(tmp_path):
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
    assert values["COUNT"] == "0"
    assert values["STATE"] == "complete_noop"
    assert Path(values["TASK_FILE"]).read_text() == ""


def test_completed_singleton_becomes_zero_task_noop(tmp_path):
    out = tmp_path / "sample_qc"
    report = out / "reports/sample_qc.tsv"
    report.parent.mkdir(parents=True)
    report.write_text("sample\tqc_status\nS1\tPASS\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "sample_variant_filtering:\n"
        "  enabled: true\n"
        f"  output_dir: {out}\n"
    )

    result = subprocess.run(
        [sys.executable, str(MANIFEST), "sample_variant_filtering", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    values = parse_values(result.stdout)
    assert values["COUNT"] == "0"
    assert values["STATE"] == "complete_noop"

    wrapper = subprocess.run(
        ["bash", str(WRAPPER), "--dry-run-submit", "sample_variant_filtering", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert wrapper.returncode == 0, wrapper.stderr
    assert "DRY RUN:" not in wrapper.stdout
    assert "state=complete_noop" in wrapper.stderr


def write_runtime_config(tmp_path: Path, mapping_rows: str) -> Path:
    mapping = tmp_path / "codon_map.tsv"
    mapping.write_text("sample\treference_key\n" + mapping_rows)
    inputs = tmp_path / "lifted"
    inputs.mkdir(exist_ok=True)
    (inputs / "S1.lifted.raw.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    config = tmp_path / "runtime.yaml"
    config.write_text(
        "codon_match:\n"
        "  paths:\n"
        f"    sample_reference_map: {mapping}\n"
        f"    input_vcf_dir: {inputs}\n"
        "  settings:\n"
        "    input_vcf_pattern: \"{sample}.lifted.raw.vcf\"\n"
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
        env={**__import__("os").environ, "SLURM_ARRAY_TASK_ID": "1"},
    )
    assert result.returncode == 0, result.stderr
    assert "runtime_status=skipped" in result.stderr
    assert "sample_not_in_pass_production_codon_map" in result.stderr
