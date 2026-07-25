import pytest

from qc_analysis.scripts.run_codon_match import (
    annotate, determine_status, is_resolved_codon, is_valid_iupac_codon, load_codon_index,
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
