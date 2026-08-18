from pathlib import Path

from qc_analysis.lib.mitos_rna import (
    comparable_gene,
    parse_result_mitos,
    per_base_assignments,
    reconcile_result_mitos_record,
)
from qc_analysis.scripts.build_mitos2_rrna_structure_table import (
    build_reference_rrna_structure_rows,
    build_reference_trna_structure_rows,
)


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

    assert status == "parsed_legacy_stockholm_structure"
    assert "Legacy fallback" in note
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

    assert status == "no_result_mitos_rna_structure"
    assert "SVG plots were not parsed" in note
    assert len(rows) == 10
    assert {row["struct_class"] for row in rows} == {"unknown"}


def mitos_line(seqid, feature_type, gene, start, end, strand, structure, score="1e-6"):
    return f"{seqid}\t{feature_type}\t{gene}\tmitfi\t{start}\t{end}\t{strand}\t{score}\t.\t.\t.\t{structure}\t1\n"


def test_result_mitos_parser_preserves_rna_names_and_normalizes_coordinates(tmp_path):
    path = tmp_path / "result.mitos"
    path.write_text(
        mitos_line("chrM", "rRNA", "rrnS", 647, 1600, 1, "." * 954)
        + mitos_line("chrM", "rRNA", "rrnL", 1670, 3228, 1, "." * 1559)
        + mitos_line("chrM", "tRNA", "trnL1", 40, 48, 1, "(((...)))")
        + mitos_line("chrM", "tRNA", "trnS2", 60, 68, -1, "(((...)))")
    )

    records = parse_result_mitos(path)

    assert [(r["feature_type"], r["gene_raw"], r["strand"]) for r in records] == [
        ("rRNA", "rrnS", "+"), ("rRNA", "rrnL", "+"),
        ("tRNA", "trnL1", "+"), ("tRNA", "trnS2", "-"),
    ]
    assert (records[0]["normalized_start"], records[0]["normalized_end"]) == (648, 1601)
    assert records[-1]["structure"] == "(((...)))"
    assert records[0]["source_file"] == str(path)


def test_human_like_rrna_result_mitos_intervals_are_validated_against_gff(tmp_path):
    fasta = tmp_path / "chrM.fa"
    fasta.write_text(">chrM\n" + "A" * 3300 + "\n")
    raw = tmp_path / "raw"; raw.mkdir()
    (raw / "result.mitos").write_text(
        mitos_line("chrM", "rRNA", "rrnS", 647, 1600, 1, "." * 954)
        + mitos_line("chrM", "rRNA", "rrnL", 1670, 3228, 1, "." * 1559)
    )
    (raw / "rrnS.sto").write_text("# STOCKHOLM 1.0\nchrM A\n#=GC SS_cons .\n//\n")
    ref = {"reference_key": "human-like", "coordinate_reference_fasta": str(fasta)}
    features = [
        {"feature_type": "rRNA", "gff_seqid": "chrM", "gene": "MT-RNR1", "gene_raw": "rrnS", "start": "648", "end": "1601", "strand": "+"},
        {"feature_type": "rRNA", "gff_seqid": "chrM", "gene": "MT-RNR2", "gene_raw": "rrnL", "start": "1671", "end": "3229", "strand": "+"},
    ]

    rows, status, _note = build_reference_rrna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_result_mitos_structure"
    by_gene = {gene: [row["genomic_pos"] for row in rows if row["rrna_gene"] == gene] for gene in ("MT-RNR1", "MT-RNR2")}
    assert (min(by_gene["MT-RNR1"]), max(by_gene["MT-RNR1"])) == (648, 1601)
    assert (min(by_gene["MT-RNR2"]), max(by_gene["MT-RNR2"])) == (1671, 3229)
    assert {row["structure_source"] for row in rows} == {str(raw / "result.mitos")}


def test_dot_bracket_pairing_is_reciprocal():
    assignments = per_base_assignments("(((...)))", 9, "result.mitos")
    assert [(i, assignments[i]["paired_local_pos"]) for i in (1, 2, 3, 7, 8, 9)] == [
        (1, 9), (2, 8), (3, 7), (7, 3), (8, 2), (9, 1),
    ]
    assert {assignments[i]["struct_class"] for i in (4, 5, 6)} == {"loop"}


def test_negative_strand_trna_uses_rna_five_prime_orientation(tmp_path):
    fasta = tmp_path / "chrM.fa"; fasta.write_text(">chrM\nACGTACGTAA\n")
    raw = tmp_path / "raw"; raw.mkdir()
    (raw / "result.mitos").write_text(mitos_line("chrM", "tRNA", "trnS2", 0, 8, -1, "(((...)))"))
    ref = {"reference_key": "ref", "coordinate_reference_fasta": str(fasta)}
    features = [{"feature_type": "tRNA", "gff_seqid": "chrM", "gene": "trnS2", "gene_raw": "trnS2", "start": "1", "end": "9", "strand": "-"}]

    rows, status, _note = build_reference_trna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_result_mitos_structure"
    assert [row["genomic_pos"] for row in rows] == list(range(9, 0, -1))
    assert rows[0]["paired_genomic_pos"] == 1
    assert rows[-1]["paired_genomic_pos"] == 9


def test_trna_copy_suffixes_reconcile_by_exact_interval(tmp_path):
    fasta = tmp_path / "chrM.fa"; fasta.write_text(">chrM\n" + "A" * 17100 + "\n")
    raw = tmp_path / "raw"; raw.mkdir()
    (raw / "result.mitos").write_text(
        mitos_line("chrM", "tRNA", "trnF", 0, 68, 1, "." * 69)
        + mitos_line("chrM", "tRNA", "trnF", 17100, 17168, 1, "." * 69)
    )
    ref = {"reference_key": "copy-suffix", "coordinate_reference_fasta": str(fasta)}
    features = [
        {"feature_type": "tRNA", "gff_seqid": "chrM", "gene": "trnF_0", "gene_raw": "trnF_0", "start": "1", "end": "69", "strand": "+"},
        {"feature_type": "tRNA", "gff_seqid": "chrM", "gene": "trnF_1", "gene_raw": "trnF_1", "start": "17101", "end": "17169", "strand": "+"},
    ]

    rows, status, note = build_reference_trna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_result_mitos_structure", note
    assert len(rows) == 138
    assert {row["structure_source"] for row in rows} == {str(raw / "result.mitos")}
    # Both final GFF copies remain present even though their extended
    # coordinates map to the same modulo-reference genomic positions.
    assert [row["genomic_pos"] for row in rows] == [*range(1, 70), *range(1, 70)]


def test_wrapped_result_mitos_interval_uses_reference_length(tmp_path):
    fasta = tmp_path / "chrM.fa"; fasta.write_text(">chrM\n" + "A" * 16752 + "\n")
    raw = tmp_path / "raw"; raw.mkdir()
    (raw / "result.mitos").write_text(
        mitos_line("chrM", "tRNA", "trnF", 16751, 69, 1, "." * 71)
    )
    ref = {"reference_key": "origin-wrap", "coordinate_reference_fasta": str(fasta)}
    features = [{
        "feature_type": "tRNA", "gff_seqid": "chrM", "gene": "trnF", "gene_raw": "trnF",
        "start": "16752", "end": "16822", "strand": "+",
    }]

    rows, status, note = build_reference_trna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_result_mitos_structure", note
    assert [row["genomic_pos"] for row in rows] == [16752, *range(1, 71)]


def test_trna_isoacceptor_numbers_remain_distinct():
    assert comparable_gene("trnL1", "tRNA") != comparable_gene("trnL2", "tRNA")
    assert comparable_gene("trnS1", "tRNA") != comparable_gene("trnS2", "tRNA")


def test_result_mitos_interval_mismatch_is_not_silently_shifted(tmp_path):
    path = tmp_path / "result.mitos"
    path.write_text(mitos_line("chrM", "rRNA", "rrnS", 647, 1600, 1, "." * 954))
    record = parse_result_mitos(path)[0]
    feature = {"feature_type": "rRNA", "gene_raw": "rrnS", "start": "649", "end": "1602", "strand": "+"}

    match, status, note = reconcile_result_mitos_record(record, [feature])

    assert match is None
    assert status == "result_mitos_gff_interval_mismatch"
    assert "648-1601" in note and "649-1602" in note


def test_circular_rna_interval_maps_across_origin(tmp_path):
    fasta = tmp_path / "chrM.fa"; fasta.write_text(">chrM\nACGTACGTAA\n")
    raw = tmp_path / "raw"; raw.mkdir()
    (raw / "result.mitos").write_text(mitos_line("chrM", "rRNA", "rrnS", 7, 11, 1, "((.))"))
    ref = {"reference_key": "circular", "coordinate_reference_fasta": str(fasta)}
    features = [{"feature_type": "rRNA", "gff_seqid": "chrM", "gene": "MT-RNR1", "gene_raw": "rrnS", "start": "8", "end": "12", "strand": "+"}]

    rows, status, _note = build_reference_rrna_structure_rows(ref, features, fasta, raw)

    assert status == "parsed_result_mitos_structure"
    assert [row["genomic_pos"] for row in rows] == [8, 9, 10, 1, 2]
    assert rows[0]["paired_genomic_pos"] == 2


def test_structure_builder_reports_length_and_gff_mismatch_statuses(tmp_path):
    fasta = tmp_path / "chrM.fa"; fasta.write_text(">chrM\nACGTACGTAA\n")
    ref = {"reference_key": "ref", "coordinate_reference_fasta": str(fasta)}
    feature = {"feature_type": "rRNA", "gff_seqid": "chrM", "gene": "MT-RNR1", "gene_raw": "rrnS", "start": "1", "end": "9", "strand": "+"}

    length_raw = tmp_path / "length"; length_raw.mkdir()
    (length_raw / "result.mitos").write_text(mitos_line("chrM", "rRNA", "rrnS", 0, 8, 1, "((..))"))
    rows, status, _note = build_reference_rrna_structure_rows(ref, [feature], fasta, length_raw)
    assert status == "result_mitos_structure_length_mismatch"
    assert {row["struct_class"] for row in rows} == {"unknown"}

    interval_raw = tmp_path / "interval"; interval_raw.mkdir()
    (interval_raw / "result.mitos").write_text(mitos_line("chrM", "rRNA", "rrnS", 1, 9, 1, "(((...)))"))
    _rows, status, note = build_reference_rrna_structure_rows(ref, [feature], fasta, interval_raw)
    assert status == "result_mitos_gff_interval_mismatch"
    assert "normalized result.mitos=2-10" in note
