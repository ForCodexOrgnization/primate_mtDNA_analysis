import pytest

from qc_analysis.scripts.run_codon_match import (
    annotate, determine_status, is_resolved_base, is_resolved_codon, is_valid_iupac_base,
    is_valid_iupac_codon, load_codon_index,
)


def row(gene='MT-CO1', codon='GCT', phase='1', strand='+'):
    return {'gene': gene, 'codon_seq': codon, 'codon_pos_in_triplet': phase,
            'strand': strand, 'ref_base_genome': 'G'}


@pytest.mark.parametrize('codon', ['GNT', 'NNN', 'AAR', 'YTG', 'gnt'])
def test_valid_iupac_codons(codon):
    assert is_valid_iupac_codon(codon)


@pytest.mark.parametrize('codon', ['', 'AT', 'ATGC', 'A-X', 'A?G'])
def test_malformed_codons(codon):
    assert not is_valid_iupac_codon(codon)


def test_production_gnt_loads_and_is_preserved(tmp_path):
    table = tmp_path / 'source.tsv'
    table.write_text('reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\tannotation_source\n'
                     'ref\t219086\tMT-CO1\t+\t1\tgnt\tG\tGenBank\n')
    index = load_codon_index(table, 'reference_key')
    assert index[('ref', 219086)][0].codon_seq == 'GNT'
    assert len(index.ambiguous_details) == 1


@pytest.mark.parametrize('base', ['N', 'R', 'Y', 'W', 'n'])
def test_valid_iupac_reference_bases(base):
    assert is_valid_iupac_base(base)
    assert not is_resolved_base(base)


def test_production_ambiguous_reference_regression(tmp_path):
    table = tmp_path / 'source.tsv'
    table.write_text('reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
                     'KY117600.1\t6903\tMT-CO1\t+\t2\tGNT\tn\n')
    index = load_codon_index(table, 'reference_key')
    annotation = index[('KY117600.1', 6903)][0]
    assert annotation.ref_base_genome == 'N'
    assert len(index.ambiguous_ref_details) == 1


def test_invalid_reference_base_remains_fatal(tmp_path):
    table = tmp_path / 'source.tsv'
    table.write_text('reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
                     'ref\t1\tMT-CO1\t+\t1\tGCT\tX\n')
    with pytest.raises(SystemExit, match='valid IUPAC'):
        load_codon_index(table, 'reference_key')


def test_only_conflicting_resolved_reference_bases_are_fatal(tmp_path):
    header = 'reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
    mixed = tmp_path / 'mixed.tsv'
    mixed.write_text(header + 'ref\t1\tMT-CO1\t+\t1\tGCT\tG\nref\t1\tMT-ND1\t+\t1\tGNT\tN\n')
    assert len(load_codon_index(mixed, 'reference_key')[('ref', 1)]) == 2
    conflict = tmp_path / 'conflict.tsv'
    conflict.write_text(header + 'ref\t1\tMT-CO1\t+\t1\tGCT\tA\nref\t1\tMT-ND1\t+\t1\tGCT\tG\n')
    with pytest.raises(SystemExit, match='Inconsistent ref_base_genome'):
        load_codon_index(conflict, 'reference_key')


def test_invalid_iupac_row_remains_fatal(tmp_path):
    table = tmp_path / 'source.tsv'
    table.write_text('reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
                     'ref\t1\tMT-CO1\t+\t1\tG-X\tG\n')
    with pytest.raises(SystemExit, match='valid IUPAC'):
        load_codon_index(table, 'reference_key')


def test_ambiguous_source_has_no_alternate_or_match():
    values, _ = annotate([row(codon='GNT')], [row(codon='GCT')], 'A', True, 'G')
    assert values['MTCODON_PRIMATE_ALT_CODON'] == '.'
    assert values['MTCODON_SOURCE_CODON_RESOLVED'] == 'no'
    assert values['MTCODON_MATCH'] == 'no'


def test_ambiguous_source_reference_is_unknown_not_mismatch():
    source = row(); source['ref_base_genome'] = 'N'
    values, candidates = annotate([source], [row()], 'A', True, 'G')
    assert values['MTCODON_SOURCE_REF_MATCH'] == 'NA'
    assert values['MTCODON_SOURCE_REF_RESOLVED'] == 'no'
    assert values['MTCODON_PRIMATE_ALT_CODON'] == '.'
    assert determine_status(1, 1, [source], [row()], values, candidates) == 'AMBIGUOUS_SOURCE_REF'


def test_resolved_source_reference_mismatch_and_match():
    mismatch = row(); mismatch['ref_base_genome'] = 'A'
    values, candidates = annotate([mismatch], [row()], 'A', True, 'G')
    assert determine_status(1, 1, [mismatch], [row()], values, candidates) == 'SOURCE_REF_MISMATCH'
    matching = row(); matching['ref_base_genome'] = 'G'
    values, candidates = annotate([matching], [row()], 'A', True, 'G')
    assert values['MTCODON_SOURCE_REF_MATCH'] == 'yes'
    assert determine_status(1, 1, [matching], [row()], values, candidates) == 'PASS'


def test_resolved_overlap_passes_despite_ambiguous_reference_candidate():
    ambiguous = row('MT-CO1'); ambiguous['ref_base_genome'] = 'N'
    resolved = row('MT-ND1'); resolved['ref_base_genome'] = 'G'
    human = [row('MT-CO1'), row('MT-ND1')]
    values, candidates = annotate([ambiguous, resolved], human, 'A', True, 'G')
    assert values['MTCODON_MATCH'] == 'yes'
    assert values['MTCODON_PRIMATE_GENE'] == 'MT-ND1'
    assert values['MTCODON_ANY_RESOLVED_SOURCE_REF'] == 'yes'
    assert determine_status(1, 1, [ambiguous, resolved], human, values, candidates) == 'PASS'


def test_ambiguous_human_has_no_match():
    values, _ = annotate([row(codon='GCT')], [row(codon='AAR')], 'A', True, 'G')
    assert values['MTCODON_HUMAN_CODON_RESOLVED'] == 'no'
    assert values['MTCODON_MATCH'] == 'no'


def test_resolved_overlap_can_pass_despite_ambiguous_annotation():
    source = [row('MT-CO1', 'GNT'), row('MT-ND1', 'GCT')]
    human = [row('MT-CO1', 'GCT'), row('MT-ND1', 'GCT')]
    values, _ = annotate(source, human, 'A', True, 'G')
    assert values['MTCODON_MATCH'] == 'yes'
    assert values['MTCODON_PRIMATE_GENE'] == 'MT-ND1'
    assert values['MTCODON_ANY_RESOLVED_PAIR'] == 'yes'


def test_resolved_helper_rejects_valid_ambiguity():
    assert is_resolved_codon('GCT')
    assert not is_resolved_codon('GNT')


@pytest.mark.parametrize('source,human,expected', [
    (row(codon='GNT'), row(codon='GCT'), 'AMBIGUOUS_CODON'),
    (row('MT-CO1', 'GNT'), row('MT-ND1', 'GNT'), 'GENE_MISMATCH'),
    (row(codon='GNT', phase='1'), row(codon='GNT', phase='2'), 'PHASE_MISMATCH'),
])
def test_ambiguity_status_does_not_mask_gene_or_phase(source, human, expected):
    values, candidates = annotate([source], [human], 'A', True, 'G')
    assert determine_status(1, 1, [source], [human], values, candidates) == expected
