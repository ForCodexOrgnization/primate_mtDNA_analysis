import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from qc_analysis.scripts.run_trna_match import (
    index, normalize_chrom, oriented, resolve_coordinate_reference_fasta,
    sample_reference_key,
)

ROOT = Path(__file__).resolve().parents[2]
HEADER = ('chrom\tpos\ttrna_id\tlocal_pos\tstruct_class\tstruct_element\tpair_type\t'
          'pair_state\tpaired_local_pos\tpaired_genomic_pos\tpaired_base\tstrand\tbase_orientation\n')

def row(chrom='species', pos=10, ident='S', strand='+', paired='G'):
    return f'{chrom}\t{pos}\t{ident}\t1\tstem\tacceptor\tWC\tpaired\t2\t20\t{paired}\t{strand}\tgenomic\n'

def fixture(tmp_path, species_rows='', human_rows='', settings='', compressed=False):
    human=tmp_path/'human.tsv'; species=tmp_path/'S1.tsv'; human.write_text(HEADER+human_rows); species.write_text(HEADER+species_rows)
    vcf='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrM\t10\t.\tA\tT\t.\tPASS\tSRC_CHROM=species;SRC_POS=10;SRC_ALT=T\n'
    inp=tmp_path/('input.vcf.gz' if compressed else 'input.vcf')
    (gzip.open(inp,'wt') if compressed else inp.open('w')).write(vcf)
    config=tmp_path/'config.yaml'; config.write_text(f'''trna_match:
  paths:
    input_vcf_dir: {tmp_path}
    fallback_input_vcf_dir: {tmp_path}
    output_dir: {tmp_path}
    reports_dir: {tmp_path}/reports
    coordinate_map_dir: {tmp_path}/maps
    human_trna_index: {human}
    species_trna_index_dir: {tmp_path}
    species_trna_index_template: "{{species_trna_index_dir}}/{{sample}}.tsv"
  settings:
    input_vcf_pattern: "S1.vcf"
    fallback_input_vcf_pattern: "S1.vcf"
    output_suffix: .out.vcf
    species_trna_lookup_ignore_chrom: false
    human_trna_lookup_ignore_chrom: false
{settings}''')
    return config,inp,tmp_path/'out.vcf'

def run(config, inp, out):
    return subprocess.run([sys.executable,str(ROOT/'qc_analysis/scripts/run_trna_match.py'),'--config',str(config),'--sample','S1','--input',str(inp),'--output',str(out)],cwd=ROOT,text=True,capture_output=True)

def status(out):
    line=next(x for x in out.read_text().splitlines() if not x.startswith('#'))
    return dict(x.split('=',1) for x in line.split('\t')[7].split(';') if '=' in x)['MTTRNA_STATUS']

@pytest.mark.parametrize('species_rows,human_rows,expected',[
    ('',row('chrM',ident='H'),'NO_SPECIES_TRNA'),
    (row(),'', 'NO_HUMAN_TRNA'),
    ('','', 'NO_SPECIES_OR_HUMAN_TRNA'),
    (row(),row('chrM',ident='H'),'OK')])
def test_missing_annotations_write_all_statuses(tmp_path,species_rows,human_rows,expected):
    config,inp,out=fixture(tmp_path,species_rows,human_rows); result=run(config,inp,out)
    assert result.returncode==0,result.stderr; assert status(out)==expected

def test_chromosome_mismatch_and_controlled_fallback(tmp_path):
    path=tmp_path/'i.tsv'; path.write_text(HEADER+row('MT'))
    idx=index(path)
    assert idx.lookup('chrM',10,False)==(None,'chromosome_mismatch')
    assert idx.lookup('chrM',10,True)[1]=='position'
    assert normalize_chrom('chrM','mitochondrial_alias')=='MT'

def test_ambiguous_position_and_duplicate_validation(tmp_path):
    path=tmp_path/'i.tsv'; path.write_text(HEADER+row('a')+row('b'))
    idx=index(path); assert idx.lookup('x',10,True)==(None,'ambiguous')
    path.write_text(HEADER+row('a')+row('a'))
    with pytest.warns(RuntimeWarning): assert index(path).duplicate_keys=={('a',10)}
    path.write_text(HEADER+row('a')+row('a',ident='conflict'))
    idx=index(path)
    assert idx.lookup('a',10)==(None,'ambiguous')
    assert idx.n_overlapping_positions==idx.n_multi_trna_positions==1

def test_negative_strand_genomic_gu_wobble_orientation():
    record={'strand':'-'}
    # genomic C/A complements to transcript G/U
    from qc_analysis.lib.match_utils import pair_type
    assert pair_type(oriented('C',record),oriented('A',record))=='GU_wobble'

def test_compressed_input_and_src_alt_flip(tmp_path):
    config,inp,out=fixture(tmp_path,row(paired='A'),row('chrM',ident='H',paired='A'),compressed=True)
    with gzip.open(inp,'rt') as handle: content=handle.read()
    with gzip.open(inp,'wt') as handle: handle.write(content.replace('\tA\tT\t.\tPASS', '\tA\tC\t.\tPASS'))
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    assert 'MTTRNA_S_ALT_PAIR_TYPE=WC' in out.read_text()  # SRC_ALT=T despite lifted ALT=C.
    assert 'MTTRNA_H_ALT_PAIR_TYPE=non_WC' in out.read_text()

def test_missing_coordinate_map_is_counted(tmp_path):
    config,inp,out=fixture(tmp_path,row(),row('chrM',ident='H'))
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    summary=(tmp_path/'reports/S1.trna_match_summary.tsv').read_text().splitlines()
    assert dict(zip(summary[0].split('\t'),summary[1].split('\t')))['missing_coordinate_map']=='1'

def test_missing_input_reports_both_paths(tmp_path):
    config,inp,out=fixture(tmp_path,row(),row('chrM')); inp.unlink(); result=run(config,inp,out)
    assert result.returncode!=0 and 'attempted primary' in result.stderr and 'fallback' in result.stderr

def test_pass_only_and_summary_controls(tmp_path):
    config,inp,out=fixture(tmp_path,row(),row('chrM'),settings='    pass_only: true\n    write_summary: false\n')
    inp.write_text(inp.read_text().replace('\tPASS\t','\tLowQual\t')); result=run(config,inp,out)
    assert result.returncode==0 and not any(not x.startswith('#') for x in out.read_text().splitlines())
    assert not (tmp_path/'reports').exists()

def test_v2_paired_rna_is_not_complemented_twice():
    from qc_analysis.scripts.run_trna_match import paired_rna
    assert paired_rna({'strand':'-','paired_base_rna':'U','paired_base':'A','index_format_version':'2'}) == 'U'

def test_sha_reference_key_resolves_exact_coordinate_fasta(tmp_path):
    digest='a'*64; key=f'mtref_{digest}'; exact=tmp_path/'Ref_chrM.fa'
    mapping=tmp_path/'sample_coordinate_reference_map.tsv'
    mapping.write_text('sample\treference_key\tcoordinate_reference_fasta\tcoordinate_reference_sequence_sha256\n'
                       f'S1\t{key}\t{exact}\t{digest}\n')
    assert sample_reference_key(mapping,'S1') == key
    assert resolve_coordinate_reference_fasta(mapping,'S1',key) == str(exact)

def test_sample_must_have_exactly_one_reference_key(tmp_path):
    mapping=tmp_path/'sample_coordinate_reference_map.tsv'
    mapping.write_text('sample\treference_key\tcoordinate_reference_fasta\n'
                       'S1\tmtref_'+'a'*64+'\ta.fa\n'
                       'S1\tmtref_'+'b'*64+'\tb.fa\n')
    with pytest.raises(ValueError,match='exactly one reference_key'):
        sample_reference_key(mapping,'S1')

def test_existing_valid_reference_index_needs_no_fasta_lookup(tmp_path):
    human=tmp_path/'human.tsv'; human.write_text(HEADER+row('chrM',ident='H'))
    digest='c'*64; key=f'mtref_{digest}'; indexes=tmp_path/'indexes'; indexes.mkdir()
    (indexes/f'{key}.tsv').write_text(HEADER+row())
    mapping=tmp_path/'sample_coordinate_reference_map.tsv'
    # Deliberately no coordinate_reference_fasta column: it must not be read.
    mapping.write_text(f'sample\treference_key\nS1\t{key}\n')
    inp=tmp_path/'input.vcf'; inp.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrM\t10\t.\tA\tT\t.\tPASS\tSRC_CHROM=species;SRC_POS=10\n')
    cfg=tmp_path/'config.yaml'; cfg.write_text(f'''trna_match:
  paths:
    input_vcf_dir: {tmp_path}
    fallback_input_vcf_dir: {tmp_path}
    output_dir: {tmp_path}
    reports_dir: {tmp_path}/reports
    coordinate_map_dir: {tmp_path}/maps
    sample_reference_map: {mapping}
    human_trna_index: {human}
    reference_trna_index_dir: {indexes}
    reference_trna_index_template: "{{reference_trna_index_dir}}/{{reference_key}}.tsv"
  settings:
    input_vcf_pattern: unused.vcf
    fallback_input_vcf_pattern: unused.vcf
    output_suffix: .out.vcf
    run_trnascan_if_missing: false
''')
    result=run(cfg,inp,tmp_path/'out.vcf')
    assert result.returncode == 0,result.stderr
    assert f'MTTRNA_REFERENCE_KEY={key}' in (tmp_path/'out.vcf').read_text()
