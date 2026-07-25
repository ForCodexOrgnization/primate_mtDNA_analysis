import csv
import gzip
import subprocess
import sys
from pathlib import Path

from qc_analysis.lib.match_utils import info_parse
from qc_analysis.scripts.run_codon_match import (
    annotate, find_overlapping_annotations, load_codon_index,
)


HEADER = 'reference_key\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\n'
SCRIPT = Path(__file__).parents[1] / 'scripts' / 'run_codon_match.py'


def row(gene, codon='AAA', phase='1', strand='+'):
    return {'gene': gene, 'codon_seq': codon, 'codon_pos_in_triplet': phase, 'strand': strand}


def test_streaming_plain_and_gzip_preserves_duplicate_positions(tmp_path, monkeypatch):
    content = HEADER + 'ref\t10\tMT-ATP8\t+\t1\tAAA\nref\t10\tMT-ATP6\t+\t1\tAAA\n'
    plain = tmp_path / 'codons.tsv'; plain.write_text(content)
    zipped = tmp_path / 'codons.tsv.gz'
    with gzip.open(zipped, 'wt') as handle:
        handle.write(content)
    # The codon loader must not use the materializing rows() helper.
    monkeypatch.setattr('qc_analysis.lib.match_utils.rows', lambda path: (_ for _ in ()).throw(AssertionError()))
    for path in (plain, zipped):
        index = load_codon_index(path, 'reference_key')
        assert [item['gene'] for item in index[('ref', 10)]] == ['MT-ATP8', 'MT-ATP6']


def test_atp_overlap_evaluates_all_pairs_and_pass_is_not_hidden():
    source = [row('MT-ATP8', 'CCC', '2'), row('MT-ATP6', 'AAA', '1')]
    human = [row('MT-ATP8', 'GGG', '3'), row('MT-ATP6', 'TAA', '1')]
    values, candidates = annotate(source, human, 'T', True)
    assert len(candidates) == 4
    assert values['MTCODON_N_PRIMATE_ANNOTATIONS'] == '2'
    assert values['MTCODON_N_HUMAN_ANNOTATIONS'] == '2'
    assert values['MTCODON_N_PAIR_CANDIDATES'] == '4'
    assert values['MTCODON_OVERLAPPING_CDS'] == 'yes'
    assert values['MTCODON_MATCH'] == 'yes'
    assert values['MTCODON_MATCHING_GENES'] == 'MT-ATP6'
    assert values['MTCODON_PRIMATE_GENE'] == 'MT-ATP6'


def test_nd4_overlap_and_row_order_independence():
    source = [row('MT-ND4L', 'CCC'), row('MT-ND4', 'AAA')]
    human = [row('MT-ND4L', 'CCC'), row('MT-ND4', 'AAA')]
    forward, _ = annotate(source, human, 'G', True)
    reverse, _ = annotate(list(reversed(source)), list(reversed(human)), 'G', True)
    assert forward == reverse
    assert forward['MTCODON_MATCH'] == 'yes'
    assert forward['MTCODON_MATCHING_GENES'] == 'MT-ND4,MT-ND4L'
    assert forward['MTCODON_AMBIGUOUS_BEST_MATCH'] == 'yes'
    assert forward['MTCODON_PRIMATE_GENE'] == 'MT-ND4'


def test_nonoverlap_regression_and_minus_strand_alt():
    values, candidates = annotate([row('MT-ND1', 'ACG', '2', '-')], [row('MT-ND1', 'ATG', '2')], 'A', True)
    assert len(candidates) == 1
    assert values['MTCODON_OVERLAPPING_CDS'] == 'no'
    assert values['MTCODON_PRIMATE_ALT_CODON'] == 'ATG'  # complement(A) is T
    assert values['MTCODON_SUPPORTED_SNV'] == 'yes'
    assert values['MTCODON_GENE_MATCH'] == values['MTCODON_PHASE_MATCH'] == values['MTCODON_MATCH'] == 'yes'


def test_overlap_report_has_required_columns(tmp_path):
    table = tmp_path / 'source.tsv'
    table.write_text(HEADER + 'ref\t10\tMT-ATP8\t+\t1\tAAA\nref\t10\tMT-ATP6\t+\t1\tAAA\n')
    assert find_overlapping_annotations(table, 'primate', 'reference_key') == [{
        'table': 'primate', 'reference_key': 'ref', 'position': 10,
        'genes': 'MT-ATP6,MT-ATP8', 'number_of_annotations': 2,
    }]


def test_four_record_smoke_and_summary_metrics(tmp_path):
    source_table = tmp_path / 'source.tsv'; human_table = tmp_path / 'human.tsv'
    source_table.write_text(HEADER + ''.join([
        'ref\t1\tMT-ND1\t+\t1\tAAA\n',
        'ref\t2\tMT-ATP8\t+\t1\tAAA\n', 'ref\t2\tMT-ATP6\t+\t1\tAAA\n',
        'ref\t3\tMT-ND4L\t+\t1\tCCC\n', 'ref\t3\tMT-ND4\t+\t1\tCCC\n']))
    human_table.write_text(HEADER.replace('reference_key\t', '') + ''.join([
        '1\tMT-ND1\t+\t1\tAAA\n',
        '2\tMT-ATP8\t+\t1\tAAA\n', '2\tMT-ATP6\t+\t1\tAAA\n',
        '3\tMT-ND4L\t+\t1\tCCC\n', '3\tMT-ND4\t+\t1\tCCC\n']))
    mapping = tmp_path / 'map.tsv'; mapping.write_text('sample\treference_key\nS1\tref\n')
    vcf = tmp_path / 'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(
        f'chrM\t{pos}\t.\tA\tG\t.\tPASS\tSRC_CHROM=ref;SRC_POS={pos};SRC_REF=A;SRC_ALT=G\n'
        for pos in (1, 2, 3, 4)))
    output = tmp_path / 'output.vcf'; reports = tmp_path / 'reports'; config = tmp_path / 'config.yaml'
    config.write_text(f'''codon_match:
  paths:
    reference_codon_table: {source_table}
    sample_reference_map: {mapping}
    human_codon_table: {human_table}
    reports_dir: {reports}
    output_dir: {tmp_path}
  settings:
    strict_phase_match: true
    output_suffix: .vcf
''')
    subprocess.run([sys.executable, str(SCRIPT), '--config', str(config), '--sample', 'S1',
                    '--input', str(vcf), '--output', str(output)], check=True)
    infos = [info_parse(line.split('\t')[7]) for line in output.read_text().splitlines() if not line.startswith('#')]
    assert [item['MTCODON_STATUS'] for item in infos] == ['PASS', 'PASS', 'PASS', 'SKIPPED_NONCODING']
    with (reports / 'S1.codon_match_summary.tsv').open() as handle:
        summary = next(csv.DictReader(handle))
    assert summary['records_with_overlapping_source_cds'] == '2'
    assert summary['records_with_overlapping_human_cds'] == '2'
    assert summary['records_with_overlapping_cds'] == '2'
    assert summary['records_with_multiple_pair_candidates'] == '2'
