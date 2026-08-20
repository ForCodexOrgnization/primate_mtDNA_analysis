import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from qc_analysis.scripts.run_trna_match import (
    canonical_trna_identity, index, normalize_chrom, normalize_trna_identity, oriented,
    resolve_coordinate_reference_fasta, sample_reference_key, trna_identity_match,
)

ROOT = Path(__file__).resolve().parents[2]
HEADER = ('chrom\tpos\ttrna_id\taa\tanticodon\tlocal_pos\tstruct_class\tstruct_element\tpair_type\t'
          'pair_state\tpaired_local_pos\tpaired_genomic_pos\tpaired_base\tstrand\tbase_orientation\n')

def row(chrom='species', pos=10, ident='S', strand='+', paired='G', aa='Phe', anticodon='GAA'):
    return f'{chrom}\t{pos}\t{ident}\t{aa}\t{anticodon}\t1\tstem\tacceptor\tWC\tpaired\t2\t20\t{paired}\t{strand}\tgenomic\n'

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

def annotations(out):
    line=next(x for x in out.read_text().splitlines() if not x.startswith('#'))
    return dict(x.split('=',1) for x in line.split('\t')[7].split(';') if '=' in x)

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

def test_ambiguous_stem_does_not_gain_deterministic_effect_or_strict_match(tmp_path):
    ambiguous=row().replace('\tWC\tpaired\t','\tambiguous\tNA\t').replace('\tG\t+\t','\tR\t+\t')
    config,inp,out=fixture(tmp_path,ambiguous,ambiguous.replace('species','chrM').replace('\tS\t','\tH\t'))
    maps=tmp_path/'maps'; maps.mkdir(); (maps/'S1.coordinate_map.tsv').write_text(
        'species_pos_original\thuman_pos_canonical\n20\t20\n')
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    text=out.read_text()
    assert 'MTTRNA_S_ALT_EFFECT=NA' in text
    assert 'MTTRNA_ALLELE_EFFECT_MATCH=.' in text
    assert 'MTTRNA_COMPENSATED=.' in text
    assert 'MTTRNA_STRICT_MATCH=no' in text

def test_unknown_identity_prevents_loop_strict_match(tmp_path):
    loop=row(aa='Leu',anticodon='NNN').replace('\tstem\tacceptor\tWC\tpaired\t2\t20\tG\t',
                       '\tloop\tanticodon_loop\tNA\tNA\t.\t.\tR\t')
    config,inp,out=fixture(tmp_path,loop,loop.replace('species','chrM').replace('\tS\t','\tH\t'))
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    values=annotations(out)
    assert values['MTTRNA_ID_MATCH']=='.'
    assert values['MTTRNA_STRICT_MATCH']=='no'

@pytest.mark.parametrize('species_aa,species_anticodon,human_aa,human_anticodon,expected',[
    ('Trp','TCA','Trp','TCA','yes'), ('Trp','TCA','Ala','TGC','no'),
    ('Ser','TGA','Ser','GCT','no'), ('Leu','TAA','Leu','TAG','no'),
    ('Leu','NNN','Leu','TAA','.'),
])
def test_loop_strict_match_requires_biological_identity(
    tmp_path,species_aa,species_anticodon,human_aa,human_anticodon,expected
):
    loop=row(ident='chrM.trna7',aa=species_aa,anticodon=species_anticodon).replace('\tstem\tacceptor\tWC\tpaired\t2\t20\tG\t',
        '\tloop\tanticodon_loop\tNA\tNA\t.\t.\tR\t')
    human=row('chrM',ident='chrM.trna8',aa=human_aa,anticodon=human_anticodon).replace(
        '\tstem\tacceptor\tWC\tpaired\t2\t20\tG\t','\tloop\tanticodon_loop\tNA\tNA\t.\t.\tR\t')
    config,inp,out=fixture(tmp_path,loop,human); result=run(config,inp,out)
    assert result.returncode==0,result.stderr
    values=annotations(out); assert values['MTTRNA_ID_MATCH']==expected
    assert values['MTTRNA_STRICT_MATCH']==('yes' if expected=='yes' else 'no')

@pytest.mark.parametrize('species_aa,human_aa,expected',[
    ('Lys','Lys','yes'), ('Lys','Gly','no'),
])
def test_stem_strict_match_requires_identity(tmp_path,species_aa,human_aa,expected):
    config,inp,out=fixture(tmp_path,row(aa=species_aa),row('chrM',aa=human_aa))
    maps=tmp_path/'maps'; maps.mkdir(); (maps/'S1.coordinate_map.tsv').write_text(
        'species_pos_original\thuman_pos_canonical\n20\t20\n')
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    values=annotations(out); assert values['MTTRNA_ID_MATCH']==expected
    assert values['MTTRNA_STRICT_MATCH']==('yes' if expected=='yes' else 'no')

@pytest.mark.parametrize('require_compensated,expected',[
    (False,'yes'), (True,'no'),
])
def test_stem_strict_match_compensated_requirement_is_configurable(
    tmp_path,require_compensated,expected
):
    setting=str(require_compensated).lower()
    config,inp,out=fixture(
        tmp_path,row(aa='Lys'),row('chrM',aa='Lys'),
        settings=f'    require_compensated_for_strict_stem: {setting}\n',
    )
    inp.write_text(inp.read_text().replace(
        '\tA\tT\t.\tPASS\tSRC_CHROM=species;SRC_POS=10;SRC_ALT=T',
        '\tA\tG\t.\tPASS\tSRC_CHROM=species;SRC_POS=10;SRC_ALT=G',
    ))
    maps=tmp_path/'maps';maps.mkdir();(maps/'S1.coordinate_map.tsv').write_text(
        'species_pos_original\thuman_pos_canonical\n20\t20\n')
    result=run(config,inp,out);assert result.returncode==0,result.stderr
    values=annotations(out)
    for field in (
        'MTTRNA_ID_MATCH','MTTRNA_REGION_MATCH','MTTRNA_ELEMENT_MATCH',
        'MTTRNA_PAIR_STATUS_MATCH','MTTRNA_PAIR_POS_MATCH','MTTRNA_ALLELE_EFFECT_MATCH',
    ):
        assert values[field]=='yes'
    assert values['MTTRNA_COMPENSATED']=='no'
    assert values['MTTRNA_STRICT_MATCH']==expected

def test_identity_aliases_are_exact_and_isoacceptor_specific():
    assert normalize_trna_identity('Phe') == normalize_trna_identity('TRNF') == 'MT-TF'
    assert normalize_trna_identity('Val') == normalize_trna_identity('MT-TV') == 'MT-TV'
    assert normalize_trna_identity('Leu') is None
    assert normalize_trna_identity('Ser') is None

@pytest.mark.parametrize('aa,anticodon,expected',[
    ('Phe','GAA','MT-TF'), ('Trp','TCA','MT-TW'), ('Lys','TTT','MT-TK'),
    ('Leu','TAA','MT-TL1'), ('Leu','UAG','MT-TL2'),
    ('Ser','UGA','MT-TS1'), ('Ser','GCU','MT-TS2'),
    ('Leu','NNN',None), ('Ser','',None),
])
def test_canonical_identity_uses_aa_and_anticodon(aa,anticodon,expected):
    assert canonical_trna_identity({'aa':aa,'anticodon':anticodon}) == expected

def test_different_record_ids_can_match_same_biological_identity(tmp_path):
    species=row(ident='chrM.trna7',aa='Trp',anticodon='TCA')
    human=row('chrM',ident='chrM.trna8',aa='Trp',anticodon='TCA')
    config,inp,out=fixture(tmp_path,species,human)
    maps=tmp_path/'maps'; maps.mkdir(); (maps/'S1.coordinate_map.tsv').write_text(
        'species_pos_original\thuman_pos_canonical\n20\t20\n')
    result=run(config,inp,out); assert result.returncode==0,result.stderr
    values=annotations(out)
    assert values['MTTRNA_S_ID']=='chrM.trna7'
    assert values['MTTRNA_H_ID']=='chrM.trna8'
    assert values['MTTRNA_S_IDENTITY']==values['MTTRNA_H_IDENTITY']=='MT-TW'
    assert values['MTTRNA_ID_MATCH']=='yes'
    summary=(tmp_path/'reports/S1.trna_match_summary.tsv').read_text().splitlines()
    counts=dict(zip(summary[0].split('\t'),summary[1].split('\t')))
    assert counts['n_trna_id_match']=='1'
    assert counts['n_trna_id_unknown']=='0'

def test_identity_match_uses_resolved_identities():
    assert trna_identity_match('MT-TF','MT-TF') == 'yes'
    assert trna_identity_match('MT-TK','MT-TG') == 'no'
    assert trna_identity_match(None,'MT-TL1') == '.'

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
