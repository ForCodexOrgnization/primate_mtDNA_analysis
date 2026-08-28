import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "qc_analysis/scripts/run_interspecies_contamination.py"


def write_vcf(path, calls):
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for chrom, pos, ref, alt, af in calls:
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSRC_POS=999;DP=200\tDP:AF\t200:{af}\n")


def run_cohort(tmp_path, samples, calls, **settings):
    vcfs = tmp_path / "lifted"; vcfs.mkdir()
    metadata = tmp_path / "samples.tsv"
    metadata.write_text("sample\tspecies\n" + "".join(f"{s}\t{sp}\n" for s, sp in samples.items()))
    for sample, rows in calls.items():
        write_vcf(vcfs / f"{sample}.lifted.raw.vcf", rows)
    output = tmp_path / "out"
    defaults = dict(min_overlap=3, min_overlap_fraction=.5, vaf_coherence_tolerance=.03,
                    min_vaf_coherence=.7)
    defaults.update(settings)
    config = tmp_path / "qc.yaml"
    config.write_text(
        "interspecies_contamination:\n  paths:\n"
        f"    input_vcf_dir: {vcfs}\n    sample_ref_file: {metadata}\n    output_dir: {output}\n"
        "  settings:\n" + "".join(f"    {k}: {v}\n" for k, v in defaults.items()))
    result = subprocess.run([sys.executable, str(SCRIPT), "--config", str(config)],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    with (output / "reports/interspecies_contamination_report.tsv").open() as handle:
        return {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}


BASE = [("chrM", p, "A", "G", .05) for p in (101, 102, 103)]
HIGH = [(c, p, r, a, .995) for c, p, r, a, _ in BASE]


def test_clear_contamination_and_exact_human_coordinate_matching(tmp_path):
    rows = run_cohort(tmp_path, {"A":"spA", "A2":"spA", "B":"spB"},
                      {"A":BASE, "A2":[], "B":HIGH + [("chrM", 104, "A", "G", .995)]})
    assert rows["A"]["interspecies_status"] == "FAIL"
    assert rows["A"]["best_source_species"] == "spB"
    assert rows["A"]["overlap_count"] == "3"
    # SRC_POS deliberately agrees for every record; matching nevertheless uses lifted keys.
    assert rows["A"]["overlap_fraction"] == "1.000000"


def test_recipient_species_background_is_removed(tmp_path):
    rows = run_cohort(tmp_path, {"A":"spA", "A2":"spA", "B":"spB"},
                      {"A":BASE, "A2":HIGH, "B":HIGH})
    assert rows["A"]["interspecies_status"] == "PASS"
    assert rows["A"]["n_lowA_after_species_background"] == "0"


def test_singleton_recipient_species_can_only_warn(tmp_path):
    rows = run_cohort(tmp_path, {"A":"spA", "B":"spB"}, {"A":BASE, "B":HIGH})
    assert rows["A"]["interspecies_status"] == "WARN"
    assert rows["A"]["classification"] == "SINGLETON_RECIPIENT_SPECIES"


def test_vaf_incoherence_warns(tmp_path):
    low = [("chrM", p, "A", "G", af) for p, af in zip((101,102,103), (.01,.10,.20))]
    rows = run_cohort(tmp_path, {"A":"spA", "A2":"spA", "B":"spB"},
                      {"A":low, "A2":[], "B":HIGH}, min_vaf_coherence=.8)
    assert rows["A"]["classification"] == "VAF_INCOHERENT"
    assert rows["A"]["interspecies_status"] == "WARN"


def test_ambiguous_source_species_warns(tmp_path):
    rows = run_cohort(tmp_path, {"A":"spA", "A2":"spA", "B":"spB", "C":"spC"},
                      {"A":BASE, "A2":[], "B":HIGH, "C":HIGH})
    assert rows["A"]["classification"] == "AMBIGUOUS_SOURCE_SPECIES"
    assert rows["A"]["interspecies_status"] == "WARN"
