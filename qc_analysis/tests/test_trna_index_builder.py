from types import SimpleNamespace
import gzip
from qc_analysis.scripts.build_trna_position_index import build
from qc_analysis.lib.trnascan_utils import validate_trna_index
import pytest
from qc_analysis.scripts.build_all_trna_indexes import resolve_fastas, trna_chrom_normalization

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
