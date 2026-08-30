import csv
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "qc_analysis/scripts/qc_array_manifest.py"
ELIGIBILITY = ROOT / "qc_analysis/scripts/qc_sample_runtime_eligibility.py"
INTERSPECIES = ROOT / "qc_analysis/scripts/run_interspecies_contamination.py"
FINAL_FILTER = ROOT / "qc_analysis/scripts/run_final_filter.py"


def values(stdout: str) -> dict[str, str]:
    result = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def write_current_status(tmp_path: Path, sample: str = "S1") -> tuple[Path, Path]:
    collection = tmp_path / "collection"
    (collection / "reports").mkdir(parents=True, exist_ok=True)
    (collection / "reports/variant_calling_collection_summary.tsv").write_text(
        f"sample\tstatus\n{sample}\tOK\n"
    )
    liftover = tmp_path / "liftover"
    (liftover / "reports").mkdir(parents=True, exist_ok=True)
    (liftover / "reports" / f"{sample}.coordinate_liftover_qc.tsv").write_text(
        f"sample\t{sample}\nstatus\tcompleted\n"
    )
    return collection, liftover


def test_mitos2_manifest_ignores_stale_numeric_task_table(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tspecies\nS1\tSpecies_one\n")
    manifest = tmp_path / "references.tsv"
    manifest.write_text("target_species\tchrM_selection_status\nSpecies_one\tselected\n")
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    sequence = "ACGT" * 4000
    (fasta_dir / "Species_one.fa").write_text(">chrM\n" + sequence + "\n")
    stale_tasks = tmp_path / "mitos2_reference_tasks.tsv"
    stale_tasks.write_text("task_id\tstatus\n1\tcompleted\n")
    out = tmp_path / "mitos2"
    config = tmp_path / "qc.yaml"
    config.write_text(
        "mitos2_annotation:\n"
        "  paths:\n"
        f"    sample_ref_file: {samples}\n"
        f"    reference_manifest: {manifest}\n"
        f"    fasta_dir: {fasta_dir}\n"
        f"    mitos2_reference_tasks: {stale_tasks}\n"
        f"    output_dir: {out}\n"
    )
    result = subprocess.run(
        [sys.executable, str(MANIFEST), "mitos2_annotation", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    parsed = values(result.stdout)
    task = Path(parsed["TASK_FILE"]).read_text().strip()
    expected_sha = hashlib.sha256(sequence.encode()).hexdigest()
    assert task == f"reference:mtref_{expected_sha}"


def test_rrna_runtime_prepares_alias_for_raw_fallback_trna(tmp_path):
    collection, liftover = write_current_status(tmp_path)
    trna_dir = tmp_path / "trna/vcf_trna"
    trna_dir.mkdir(parents=True)
    fallback = trna_dir / "S1.lifted.trna.vcf"
    fallback.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rrna_map = tmp_path / "rrna_map.tsv"
    rrna_map.write_text("sample\treference_key\nS1\tmtref_abc\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "collect_variant_calling:\n"
        f"  outdir: {collection}\n"
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    output_dir: {liftover}\n"
        "rrna_match:\n"
        "  paths:\n"
        f"    input_vcf_dir: {trna_dir}\n"
        f"    fallback_codon_vcf_dir: {tmp_path / 'codon'}\n"
        f"    fallback_raw_vcf_dir: {liftover / 'vcf_lifted_raw'}\n"
        f"    output_dir: {tmp_path / 'rrna'}\n"
        f"    reports_dir: {tmp_path / 'rrna/reports'}\n"
        f"    sample_reference_map: {rrna_map}\n"
        "  settings:\n"
        "    input_vcf_pattern: \"{sample}.lifted.codon.trna.vcf\"\n"
        "    fallback_codon_vcf_pattern: \"{sample}.lifted.codon.vcf\"\n"
        "    fallback_raw_vcf_pattern: \"{sample}.lifted.raw.vcf\"\n"
        "    output_suffix: \".lifted.codon.trna.rrna.vcf\"\n"
    )
    result = subprocess.run(
        [sys.executable, str(ELIGIBILITY), "rrna_match", "S1", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert "ELIGIBLE=1" in result.stdout
    primary = trna_dir / "S1.lifted.codon.trna.vcf"
    assert primary.is_symlink()
    assert primary.resolve() == fallback.resolve()


def minimal_vcf(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )


def test_interspecies_ignores_stale_failed_liftover_vcf(tmp_path):
    metadata = tmp_path / "samples.tsv"
    metadata.write_text("sample\tspecies\nS1\tSpecies_one\nS2\tSpecies_two\n")
    collection = tmp_path / "collection"
    (collection / "reports").mkdir(parents=True)
    (collection / "reports/variant_calling_collection_summary.tsv").write_text(
        "sample\tstatus\nS1\tOK\nS2\tOK\n"
    )
    liftover = tmp_path / "liftover"
    (liftover / "reports").mkdir(parents=True)
    (liftover / "reports/S1.coordinate_liftover_qc.tsv").write_text("sample\tS1\nstatus\tcompleted\n")
    (liftover / "reports/S2.coordinate_liftover_qc.tsv").write_text("sample\tS2\nstatus\tskipped_anchor_validation\n")
    minimal_vcf(liftover / "vcf_lifted_raw/S1.lifted.raw.vcf")
    minimal_vcf(liftover / "vcf_lifted_raw/S2.lifted.raw.vcf")
    out = tmp_path / "interspecies"
    config = tmp_path / "qc.yaml"
    config.write_text(
        "collect_variant_calling:\n"
        f"  outdir: {collection}\n"
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    output_dir: {liftover}\n"
        "interspecies_contamination:\n"
        "  enabled: true\n"
        "  paths:\n"
        f"    input_vcf_dir: {liftover / 'vcf_lifted_raw'}\n"
        "    input_vcf_pattern: \"{sample}.lifted.raw.vcf\"\n"
        f"    sample_ref_file: {metadata}\n"
        f"    output_dir: {out}\n"
        "  settings:\n"
        "    dp_min: 100\n"
    )
    result = subprocess.run(
        [sys.executable, str(INTERSPECIES), "--config", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    report = out / "reports/interspecies_contamination_report.tsv"
    with report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample"] for row in rows] == ["S1"]
    assert "ignoring stale/non-current lifted VCFs: S2" in result.stderr


def test_final_filter_rebuild_removes_stale_final_file(tmp_path):
    collection = tmp_path / "collection"
    (collection / "reports").mkdir(parents=True)
    (collection / "reports/variant_calling_collection_summary.tsv").write_text(
        "sample\tspecies\tstatus\nS1\tSpecies_one\tOK\n"
    )
    liftover = tmp_path / "liftover"
    (liftover / "reports").mkdir(parents=True)
    (liftover / "reports/S1.coordinate_liftover_qc.tsv").write_text("sample\tS1\nstatus\tcompleted\n")
    intra = tmp_path / "intra.tsv"
    intra.write_text("sample\tcontamination_status\tqc_status\nS1\tno_strong_evidence\tPASS\n")
    sample_qc = tmp_path / "sample_qc.tsv"
    sample_qc.write_text(
        "sample\tqc_status\tfailed_criteria\tmt_median_coverage\tPercent_100\tnuclear_median_coverage\tmtcn_median\tMAD\n"
        "S1\tFAIL\tlow_mt_coverage\t10\t10\t30\t60\t0.1\n"
    )
    stale_human = tmp_path / "stale_human.tsv"
    stale_human.write_text("sample\thuman_contamination_status\nS1\tFAIL\n")
    out = tmp_path / "final"
    stale = out / "final_vcf/S1.final.vcf.gz"
    stale.parent.mkdir(parents=True)
    stale.write_text("stale\n")
    config = tmp_path / "qc.yaml"
    config.write_text(
        "coordinate_liftover:\n"
        "  paths:\n"
        f"    output_dir: {liftover}\n"
        "final_filter:\n"
        "  enabled: true\n"
        f"  collected_dir: {collection}\n"
        f"  output_dir: {out}\n"
        "  required_sample_reports: intraspecies,sample_qc\n"
        "  optional_sample_reports: interspecies\n"
        "  strict_missing_samples: true\n"
        "  sample_reports:\n"
        "    intraspecies:\n"
        f"      path: {intra}\n"
        "      status_columns: contamination_status,qc_status\n"
        "    human:\n"
        f"      path: {stale_human}\n"
        "      status_columns: human_contamination_status\n"
        "    sample_qc:\n"
        f"      path: {sample_qc}\n"
        "      status_columns: qc_status\n"
        "    interspecies:\n"
        f"      path: {tmp_path / 'missing_inter.tsv'}\n"
        "  sample_fail_status:\n"
        "    intraspecies: high_confidence_contaminated,FAIL\n"
        "    human: FAIL\n"
        "    sample_qc: FAIL\n"
        "    interspecies: FAIL\n"
        "  vcf_sources:\n"
        "    liftover:\n"
        f"      dir: {tmp_path / 'no_vcfs'}\n"
        "      pattern: \"{sample}.lifted.raw.vcf\"\n"
    )
    result = subprocess.run(
        [sys.executable, str(FINAL_FILTER), "--config", str(config)],
        cwd=ROOT, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert not stale.exists()
    report = out / "reports/final_sample_qc.tsv"
    with report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["final_sample_status"] == "FAIL"
    assert "human:" not in rows[0]["final_sample_fail_reasons"]
