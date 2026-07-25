import csv
import subprocess
import sys
from pathlib import Path

import pytest

from qc_analysis.lib.match_utils import info_parse
from qc_analysis.scripts.run_codon_match import is_supported_snv, load_sample_reference_map


SCRIPT = Path(__file__).parents[1] / 'scripts' / 'run_codon_match.py'


def make_case(tmp_path, ref='A', alt='G', *, reference_mode=True, mapped=True):
    source = tmp_path / 'source.tsv'
    key = 'reference_key' if reference_mode else 'sample'
    value = 'ref-A' if reference_mode else 'S1'
    source.write_text(f'{key}\tpos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
                      f'{value}\t10\tMT-ND1\t+\t1\tAAA\tA\n')
    human = tmp_path / 'human.tsv'
    human.write_text('pos\tgene\tstrand\tcodon_pos_in_triplet\tcodon_seq\tref_base_genome\n'
                     '10\tMT-ND1\t+\t1\tGAA\tA\n')
    vcf = tmp_path / 'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   f'chrM\t10\t.\tA\tG\t.\tPASS\tSRC_POS=10;SRC_REF={ref};SRC_ALT={alt}\n')
    reports = tmp_path / 'reports'
    paths = ''
    if reference_mode:
        mapping = tmp_path / 'map.tsv'
        mapping.write_text('sample\treference_key\n' + ('S1\tref-A\n' if mapped else 'OTHER\tref-A\n'))
        paths = f'    reference_codon_table: {source}\n    sample_reference_map: {mapping}\n'
    else:
        paths = f'    all_primate_position_codon_table: {source}\n'
    config = tmp_path / 'config.yaml'
    config.write_text(f'''codon_match:
  paths:
{paths}    human_codon_table: {human}
    reports_dir: {reports}
    output_dir: {tmp_path}
  settings:
    strict_phase_match: true
    output_suffix: .vcf
''')
    output = tmp_path / 'output.vcf'
    result = subprocess.run(
        [sys.executable, str(SCRIPT), '--config', str(config), '--sample', 'S1',
         '--input', str(vcf), '--output', str(output)], text=True, capture_output=True,
    )
    return result, output, reports


@pytest.mark.parametrize('ref,alt', [
    ('A', 'AT'), ('AT', 'A'), ('A', 'G,T'), ('A', '<DEL>'), ('A', 'N'), ('', 'G'), ('A', ''),
])
def test_non_snv_alleles_are_not_truncated(tmp_path, ref, alt):
    result, output, reports = make_case(tmp_path, ref, alt)
    assert result.returncode == 0, result.stderr
    record = next(line for line in output.read_text().splitlines() if not line.startswith('#'))
    info = info_parse(record.split('\t')[7])
    expected_status = 'SOURCE_REF_MISMATCH' if ref.strip().upper() not in {'A', 'C', 'G', 'T'} else 'UNSUPPORTED_NON_SNV'
    assert info['MTCODON_STATUS'] == expected_status
    assert info['MTCODON_SUPPORTED_SNV'] == 'no'
    assert info['MTCODON_PRIMATE_ALT_CODON'] == '.'
    assert info['MTCODON_MATCH'] == 'no'
    with (reports / 'S1.codon_match_summary.tsv').open() as handle:
        summary = next(csv.DictReader(handle, delimiter='\t'))
    assert summary[f'status_{expected_status}'] == '1'


def test_supported_snv_keeps_normal_matching(tmp_path):
    result, output, _ = make_case(tmp_path)
    assert result.returncode == 0, result.stderr
    info = info_parse(next(line for line in output.read_text().splitlines()
                           if not line.startswith('#')).split('\t')[7])
    assert info['MTCODON_STATUS'] == 'PASS'
    assert info['MTCODON_SUPPORTED_SNV'] == 'yes'
    assert info['MTCODON_PRIMATE_ALT_CODON'] == 'GAA'


def test_reference_mode_requires_sample_mapping(tmp_path):
    result, _, _ = make_case(tmp_path, mapped=False)
    assert result.returncode != 0
    assert 'missing from sample_reference_map' in result.stderr


def test_conflicting_duplicate_sample_mapping_is_fatal(tmp_path):
    mapping = tmp_path / 'map.tsv'
    mapping.write_text('sample\treference_key\nS1\tref-A\nS1\tref-B\n')
    with pytest.raises(SystemExit, match="Conflicting reference keys for sample 'S1'.*'ref-A' versus 'ref-B'"):
        load_sample_reference_map(mapping)


def test_identical_duplicate_mapping_is_tolerated(tmp_path):
    mapping = tmp_path / 'map.tsv'
    mapping.write_text('sample\treference_key\nS1\tref-A\nS1\tref-A\n')
    assert load_sample_reference_map(mapping) == {'S1': 'ref-A'}


def test_historical_sample_level_lookup_remains_supported(tmp_path):
    result, output, _ = make_case(tmp_path, reference_mode=False)
    assert result.returncode == 0, result.stderr
    assert 'MTCODON_STATUS=PASS' in output.read_text()


@pytest.mark.parametrize('ref,alt,expected', [(' a ', ' g ', True), ('A', 'G,T', False), ('N', 'G', False)])
def test_supported_snv_validation_normalizes_and_rejects(ref, alt, expected):
    assert is_supported_snv(ref, alt) is expected
