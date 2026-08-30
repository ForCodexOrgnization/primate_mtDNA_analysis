import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FINAL_FILTER = ROOT / "qc_analysis/scripts/run_final_filter.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_final_filter_streaming_test", FINAL_FILTER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reset_managed_outputs_preserves_scheduler_state(tmp_path):
    module = load_module()
    out = tmp_path / "final_filter"
    scheduler_log = out / "logs/job_arrays/123_1.err"
    scheduler_log.parent.mkdir(parents=True)
    scheduler_log.write_text("keep me\n")
    manifest = out / "job_arrays/final_filter.current.tsv"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("keep me too\n")
    stale_final = out / "final_vcf/stale.final.vcf.gz"
    stale_final.parent.mkdir(parents=True)
    stale_final.write_text("stale\n")
    stale_report = out / "reports/final_sample_qc.tsv"
    stale_report.parent.mkdir(parents=True)
    stale_report.write_text("stale\n")

    module.reset_managed_outputs(out)

    assert scheduler_log.read_text() == "keep me\n"
    assert manifest.read_text() == "keep me too\n"
    assert not stale_final.exists()
    assert not stale_report.exists()
    assert (out / "final_vcf").is_dir()
    assert (out / "reports").is_dir()


def test_variant_report_store_loads_only_requested_sample(tmp_path):
    module = load_module()
    report = tmp_path / "orthology.tsv"
    with report.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "human_chrom",
                "human_pos",
                "human_ref",
                "human_alt",
                "region_type",
                "orthology_match_status",
                "orthology_fail_reason",
            ]
        )
        writer.writerow(["S1", "chrM", "100", "A", "G", "CDS", "PASS", ""])
        writer.writerow(["S1", "chrM", "200", "C", "T", "rRNA", "FAIL", "RRNA_LOW_CONF"])
        writer.writerow(["S2", "chrM", "300", "G", "A", "tRNA", "PASS", ""])

    db = tmp_path / "evidence.sqlite3"
    connection = module.build_variant_report_store(
        db,
        {
            "orthology": {
                "path": str(report),
                "coordinate_system": "human",
                "status_columns": "orthology_match_status",
                "fail_status": "FAIL",
            }
        },
        {"S1", "S2"},
    )
    try:
        flags, annotations = module.load_sample_variant_evidence(connection, "S1")
        assert ("chrM", "300", "G", "A") not in annotations
        assert annotations[("chrM", "100", "A", "G")]["region_type"] == "CDS"
        assert annotations[("chrM", "200", "C", "T")]["orthology_fail_reason"] == "RRNA_LOW_CONF"
        assert flags[("chrM", "200", "C", "T")] == ["orthology:FAIL"]
        assert ("chrM", "100", "A", "G") not in flags
    finally:
        connection.close()
