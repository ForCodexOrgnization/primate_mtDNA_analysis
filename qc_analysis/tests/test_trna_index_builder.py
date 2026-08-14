from types import SimpleNamespace
import gzip
import csv
import hashlib
import sys
from qc_analysis.scripts.build_trna_position_index import build
from qc_analysis.lib.trnascan_utils import validate_trna_index
import pytest
from qc_analysis.scripts import build_all_trna_indexes
from qc_analysis.scripts.build_all_trna_indexes import add_coordinate_fastas, resolve_fastas, trna_chrom_normalization

def test_existing_output_index_build(tmp_path):
    fasta=tmp_path/'r.fa';fasta.write_text('>chrM\nACGT\n')
    out=tmp_path/'r.out';out.write_text('chrM 1 1 4 Phe GAA 0 0 20.0\n')
    ss=tmp_path/'r.ss';ss.write_text('chrM.trna1 (1-4)\nSeq: ACGU\nStr: ><..\n')
    dest=tmp_path/'r.tsv.gz';a=SimpleNamespace(reference_key='R',fasta=str(fasta),trnascan_out=str(out),trnascan_ss=str(ss),run_trnascan=False,output=str(dest),overwrite=False,chrom_normalization='none',max_sequence_mismatch_rate=0)
    row=build(a);assert row['n_index_rows']==4 and row['n_positive_strand_trna']==1
    with gzip.open(dest,'rt') as h:assert 'paired_base_rna' in h.readline()
    assert validate_trna_index(dest,'R')['n_rows']==4
    with pytest.raises(ValueError,match='reference_key'):
        validate_trna_index(dest,'wrong')

def test_partial_existing_index_is_rejected(tmp_path):
    corrupt=tmp_path/'partial.tsv.gz'; corrupt.write_bytes(b'\x1f\x8b\x08')
    with pytest.raises((EOFError,OSError)):
        validate_trna_index(corrupt,'R')

def test_human_and_species_normalization_differ():
    settings={'human_trna_chrom_norm':'mitochondrial_alias','species_trna_chrom_norm':'strip_chr'}
    assert trna_chrom_normalization('human',settings)=='mitochondrial_alias'
    assert trna_chrom_normalization('Pan_troglodytes',settings)=='strip_chr'

def test_fasta_resolution_rejects_wg_and_ambiguous_multicontig(tmp_path):
    mito=tmp_path/'mito.fa'; mito.write_text('>MT\nACGT\n')
    wg=tmp_path/'wg.fa'; wg.write_text('>chr1\nAAAA\n>MT\nACGT\n')
    manifest=tmp_path/'manifest.tsv'
    manifest.write_text(f'reference_key\tchrM_fasta_path\twg_expected_output_fasta\nR\t{mito}\t{wg}\n')
    assert resolve_fastas(manifest,1,10)['R']['path']==str(mito)
    manifest.write_text(f'reference_key\tfasta_path\nR\t{wg}\n')
    with pytest.raises(ValueError,match='multi-contig'):
        resolve_fastas(manifest,1,10)

def test_generic_map_keeps_reference_eligible_for_trnascan_without_codon_qc(tmp_path):
    fasta=tmp_path/'failed-mitos2.fa'; fasta.write_text('>chrM\nACGT\n')
    generic_rows=[{'sample':'S1','reference_key':'mtref_failed_codon_qc',
                   'coordinate_reference_fasta':str(fasta)}]
    resolved=add_coordinate_fastas({},generic_rows,1,10)
    assert resolved['mtref_failed_codon_qc']['path']==str(fasta)

def test_coordinate_map_accepts_identical_fasta_aliases(tmp_path):
    one=tmp_path/'one.fa'; one.write_text('>MT\nACGT\n')
    two=tmp_path/'two.fa'; two.write_text('>different_header\nac gt\n')
    sha=hashlib.sha256(b'ACGT').hexdigest(); key=f'mtref_{sha}'
    rows=[{'sample':'S1','reference_key':key,'coordinate_reference_fasta':str(one)},
          {'sample':'S2','reference_key':key,'coordinate_reference_fasta':str(two)}]
    resolved=add_coordinate_fastas({},rows,1,10)[key]
    assert resolved['path']==min(str(one),str(two))
    assert resolved['n_fasta_aliases']==2
    assert resolved['validation_status']=='VALID_IDENTICAL_FASTA_ALIASES'
    assert resolved['normalized_sequence_sha256']==sha

def test_coordinate_map_rejects_different_fasta_alias_sequences(tmp_path):
    one=tmp_path/'one.fa'; one.write_text('>MT\nACGT\n')
    two=tmp_path/'two.fa'; two.write_text('>MT\nTGCA\n')
    rows=[{'sample':'S1','reference_key':'R','coordinate_reference_fasta':str(one)},
          {'sample':'S2','reference_key':'R','coordinate_reference_fasta':str(two)}]
    with pytest.raises(ValueError,match='FAIL_HASH_CONFLICT'):
        add_coordinate_fastas({},rows,1,10)

def test_coordinate_map_rejects_reference_key_hash_mismatch(tmp_path):
    fasta=tmp_path/'one.fa'; fasta.write_text('>MT\nACGT\n')
    rows=[{'sample':'S1','reference_key':f'mtref_{"0"*64}',
           'coordinate_reference_fasta':str(fasta)}]
    with pytest.raises(ValueError,match='FAIL_REFERENCE_KEY_HASH_MISMATCH'):
        add_coordinate_fastas({},rows,1,10)

def test_coordinate_map_deduplicates_references_and_separates_human(tmp_path,monkeypatch):
    one=tmp_path/'one.fa';one.write_text('>MT\nACGT\n')
    one_alias=tmp_path/'one-alias.fa';one_alias.write_text('>another_name\nACGT\n')
    two=tmp_path/'two.fa';two.write_text('>MT\nTGCA\n')
    sample_map=tmp_path/'map.tsv'
    sample_map.write_text(
        'sample\treference_key\tcoordinate_reference_fasta\tcoordinate_reference_sequence_sha256\n'
        f'S1\tR1\t{one}\t{hashlib.sha256(b"ACGT").hexdigest()}\n'
        f'S2\tR1\t{one_alias}\t{hashlib.sha256(b"ACGT").hexdigest()}\n'
        f'S3\tR2\t{two}\t{hashlib.sha256(b"TGCA").hexdigest()}\n')
    manifest=tmp_path/'tasks.tsv';config=tmp_path/'config.yaml'
    config.write_text(
        'trna_match:\n  paths:\n'
        f'    sample_reference_map: {sample_map}\n    human_fasta: {tmp_path}/human.fa\n'
        f'    human_trna_index: {tmp_path}/human.tsv.gz\n    reference_trna_index_dir: {tmp_path}/indexes\n'
        '    reference_trna_index_template: "{reference_trna_index_dir}/{reference_key}.tsv.gz"\n'
        f'    index_build_reports_dir: {tmp_path}/reports\n    trnascan_output_dir: {tmp_path}/scan\n'
        '  settings:\n    min_mitochondrial_reference_length: 1\n    max_mitochondrial_reference_length: 10\n')
    monkeypatch.setattr(sys,'argv',['build_all_trna_indexes.py','--config',str(config),'--task-manifest',str(manifest)])
    build_all_trna_indexes.main()
    with manifest.open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
    assert [row['reference_key'] for row in rows]==['human','R1','R2']
    assert rows[0]['fasta']==str(tmp_path/'human.fa')
    assert rows[1]['n_fasta_aliases']=='2'
    assert rows[1]['validation_status']=='VALID_IDENTICAL_FASTA_ALIASES'

@pytest.mark.parametrize('rows,message',[
    ([('S1','R','', '')],'missing coordinate_reference_fasta'),
])
def test_coordinate_map_rejects_missing_fasta(tmp_path,rows,message):
    for name in ('one.fa','two.fa'):(tmp_path/name).write_text('>MT\nACGT\n')
    mapped=[{'sample':sample,'reference_key':key,'coordinate_reference_fasta':str(tmp_path/path) if path else '',
             'coordinate_reference_sequence_sha256':sha} for sample,key,path,sha in rows]
    with pytest.raises(ValueError,match=message):add_coordinate_fastas({},mapped,1,10)

def test_coordinate_map_rejects_sequence_hash_mismatch(tmp_path):
    fasta=tmp_path/'one.fa';fasta.write_text('>MT\nACGT\n')
    rows=[{'sample':'S1','reference_key':'R','coordinate_reference_fasta':str(fasta),
           'coordinate_reference_sequence_sha256':'0'*64}]
    with pytest.raises(ValueError,match='SHA256 mismatch'):add_coordinate_fastas({},rows,1,10)
