from pathlib import Path

from qc_analysis.scripts.build_mitos2_rrna_structure_table import build_reference_rrna_structure_rows


def reference(tmp_path):
    fasta = tmp_path / "chrM.fa"
    fasta.write_text(">chrM\nGAAAAAAAAC\n")
    raw = tmp_path / "raw"
    raw.mkdir()
    ref = {
        "reference_key": "ref1",
        "reference_species": "Species one",
        "coordinate_reference_accession": "ACC.1",
        "coordinate_reference_fasta": str(fasta),
        "coordinate_reference_sequence_sha256": "sha1",
    }
    features = [{
        "feature_type": "rRNA",
        "gff_seqid": "chrM",
        "gene": "MT-RNR1",
        "gene_raw": "rrnS",
        "start": "1",
        "end": "10",
        "strand": "+",
    }]
    return ref, fasta, raw, features


def test_stockholm_ss_cons_builds_explicit_pair_rows(tmp_path):
    ref, fasta, raw, features = reference(tmp_path)
    (raw / "rrnS.sto").write_text(
        "# STOCKHOLM 1.0\n"
        "chrM GAAAAAAAAC\n"
        "#=GC SS_cons <........>\n"
        "//\n"
    )

    rows, status, note = build_reference_rrna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_machine_readable_structure"
    assert "Stockholm/Infernal" in note
    assert len(rows) == 10
    assert rows[0]["struct_class"] == "stem"
    assert rows[0]["paired_genomic_pos"] == 10
    assert rows[0]["paired_local_pos"] == 10
    assert rows[0]["paired_base"] == "C"
    assert rows[0]["pair_type"] == "G-C"
    assert rows[0]["pair_state"] == "canonical"
    assert rows[1]["struct_class"] == "loop"
    assert rows[1]["pair_state"] == "unpaired"
    assert rows[-1]["paired_genomic_pos"] == 1


def test_missing_machine_readable_structure_keeps_unknown_rows(tmp_path):
    ref, fasta, raw, features = reference(tmp_path)
    (raw / "rrnS.svg").write_text("<svg></svg>\n")

    rows, status, note = build_reference_rrna_structure_rows(ref, features, fasta, raw)

    assert status == "no_machine_readable_structure"
    assert "SVG plots were not parsed" in note
    assert len(rows) == 10
    assert {row["struct_class"] for row in rows} == {"unknown"}
