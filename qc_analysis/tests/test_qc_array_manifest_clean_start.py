#!/usr/bin/env python3
import csv
import hashlib
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "qc_analysis" / "scripts" / "qc_array_manifest.py"


def run_manifest(step: str, config: Path) -> tuple[Path, Path, str]:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), step, str(config)],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    values = {}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return Path(values["TASK_FILE"]), Path(values["MANIFEST"]), completed.stderr


def write_minimal_config(tmp_path: Path) -> Path:
    sample_ref = tmp_path / "sample_ref.tsv"
    sample_ref.write_text("sample\tspecies\nS1\tSpecies_one\n")

    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    fasta = fasta_dir / "Species_one.fa"
    fasta.write_text(">chrM\nACGTACGT\n")

    reference_manifest = tmp_path / "reference_manifest.tsv"
    reference_manifest.write_text(
        "target_species\tfinal_chrM_species\tchrM_expected_output_fasta\tchrM_selection_status\n"
        f"Species_one\tSpecies_one\t{fasta}\tselected\n"
    )

    config = tmp_path / "qc.yaml"
    config.write_text(
        f"""
coordinate_liftover:
  paths:
    sample_ref_file: {sample_ref}
    output_dir: {tmp_path / 'liftover'}
mitos2_annotation:
  paths:
    sample_ref_file: {sample_ref}
    reference_manifest: {reference_manifest}
    fasta_dir: {fasta_dir}
    output_dir: {tmp_path / 'mitos2'}
    mitos2_reference_tasks: {tmp_path / 'mitos2' / 'mitos2_reference_tasks.tsv'}
    sample_coordinate_reference_map: {tmp_path / 'sample_coordinate_reference_map.tsv'}
codon_match:
  paths:
    input_vcf_dir: {tmp_path / 'liftover' / 'vcf_lifted_raw'}
    output_dir: {tmp_path / 'codon'}
    sample_reference_map: {tmp_path / 'mitos2' / 'codon_sample_reference_map.tsv'}
  settings:
    input_vcf_pattern: "{{sample}}.lifted.raw.vcf"
    output_suffix: ".lifted.codon.vcf"
trna_match:
  paths:
    input_vcf_dir: {tmp_path / 'codon' / 'vcf_codon'}
    fallback_input_vcf_dir: {tmp_path / 'liftover' / 'vcf_lifted_raw'}
    output_dir: {tmp_path / 'trna'}
    sample_reference_map: {tmp_path / 'sample_coordinate_reference_map.tsv'}
  settings:
    input_vcf_pattern: "{{sample}}.lifted.codon.vcf"
    fallback_input_vcf_pattern: "{{sample}}.lifted.raw.vcf"
    output_suffix: ".lifted.codon.trna.vcf"
rrna_match:
  paths:
    input_vcf_dir: {tmp_path / 'trna' / 'vcf_trna'}
    fallback_codon_vcf_dir: {tmp_path / 'codon' / 'vcf_codon'}
    fallback_raw_vcf_dir: {tmp_path / 'liftover' / 'vcf_lifted_raw'}
    output_dir: {tmp_path / 'rrna'}
    sample_reference_map: {tmp_path / 'mitos2' / 'sample_coordinate_reference_map.tsv'}
  settings:
    input_vcf_pattern: "{{sample}}.lifted.codon.trna.vcf"
    fallback_codon_vcf_pattern: "{{sample}}.lifted.codon.vcf"
    fallback_raw_vcf_pattern: "{{sample}}.lifted.raw.vcf"
    output_suffix: ".lifted.codon.trna.rrna.vcf"
""".strip()
        + "\n"
    )
    return config


def test_mitos2_clean_start_uses_reference_key_not_sample_id(tmp_path):
    config = write_minimal_config(tmp_path)
    task_file, manifest_file, _ = run_manifest("mitos2_annotation", config)

    sequence_sha = hashlib.sha256(b"ACGTACGT").hexdigest()
    assert task_file.read_text().splitlines() == [f"reference:mtref_{sequence_sha}"]

    with manifest_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["item_type"] == "reference"
    assert rows[0]["item"].startswith("reference:mtref_")
    assert rows[0]["item"] != "reference:S1"


def test_codon_clean_start_keeps_future_input_for_runtime_revalidation(tmp_path):
    config = write_minimal_config(tmp_path)
    task_file, manifest_file, stderr = run_manifest("codon_match", config)

    assert task_file.read_text().splitlines() == ["S1"]
    with manifest_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["status_at_submission"] == "runtime_revalidate"
    assert rows[0]["expected_input"].endswith("S1.lifted.raw.vcf")
    assert "missing inputs 1" in stderr
