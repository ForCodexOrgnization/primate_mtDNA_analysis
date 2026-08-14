import csv
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/run_mitos2_annotation.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_mitos2_mapping", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_references(tmp_path, sequence_b="TTTT", origin_b="Species_A"):
    module = load_module()
    fasta_dir = tmp_path / "Ref_chrM"
    fasta_dir.mkdir()
    (fasta_dir / "Species_A.fa").write_text(">A\nACGT\n")
    (fasta_dir / "Species_B.fa").write_text(f">B\n{sequence_b}\n")
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, [
        {"target_species": "Species_A", "final_chrM_species": "Species_A", "final_chrM_accession": "ACC_A"},
        {"target_species": "Species_B", "final_chrM_species": origin_b, "final_chrM_accession": "ACC_B"},
    ])
    samples = tmp_path / "samples.tsv"
    write_tsv(samples, [
        {"sample": "sample_a", "species": "Species_A"},
        {"sample": "sample_b", "species": "Species_B"},
    ])
    paths = {
        "reference_manifest": str(manifest),
        "sample_ref_file": str(samples),
        "final_chrM_fasta_dir": str(fasta_dir),
        "mitos2_raw_dir": str(tmp_path / "raw"),
    }
    return module, module.references(paths)


def test_reference_provenance_does_not_cross_link_samples(tmp_path):
    module, refs = build_references(tmp_path)
    rows = module.sample_reference_rows(refs)

    assert {row["sample"] for row in rows} == {"sample_a", "sample_b"}
    assert len(rows) == 2
    by_sample = {row["sample"]: row for row in rows}
    assert by_sample["sample_a"]["coordinate_reference_fasta"].endswith("Species_A.fa")
    assert by_sample["sample_b"]["coordinate_reference_fasta"].endswith("Species_B.fa")
    assert by_sample["sample_a"]["reference_key"] != by_sample["sample_b"]["reference_key"]
    # Species_B's FASTA biologically originated in Species_A, but sample_a is
    # still linked exclusively to Species_A's actual variant-calling FASTA.
    ref_b = next(ref for ref, _ in refs if ref["target_species"] == "Species_B")
    assert ref_b["reference_species"] == "Species_A"


def test_exact_sequences_share_hash_identity_without_duplicate_samples(tmp_path):
    module, refs = build_references(tmp_path, sequence_b="ACGT")
    rows = module.sample_reference_rows(refs)

    assert len(refs) == 1
    assert len({row["reference_key"] for row in rows}) == 1
    assert len({row["sample"] for row in rows}) == 2
    assert rows[0]["reference_key"].startswith("mtref_")


def test_mapping_validation_rejects_distinct_keys_with_full_diagnostic():
    module = load_module()
    rows = [
        {"sample": "ERS1", "species": "Species_A", "reference_key": "mtref_a",
         "coordinate_reference_fasta": "A.fa", "coordinate_reference_accession": "A",
         "coordinate_reference_sequence_sha256": "a"},
        {"sample": "ERS1", "species": "Species_A", "reference_key": "mtref_b",
         "coordinate_reference_fasta": "B.fa", "coordinate_reference_accession": "B",
         "coordinate_reference_sequence_sha256": "b"},
    ]

    with pytest.raises(SystemExit) as error:
        module.validate_sample_reference_rows(rows, "codon_sample_reference_map.tsv")
    message = str(error.value)
    assert all(value in message for value in ("ERS1", "Species_A", "mtref_a", "mtref_b", "A.fa", "B.fa", "SHA256"))


def test_mapping_validation_collapses_exact_duplicates():
    module = load_module()
    row = {field: "value" for field in module.SAMPLE_REFERENCE_FIELDS}
    assert module.validate_sample_reference_rows([row, dict(row)], "map.tsv") == [row]
