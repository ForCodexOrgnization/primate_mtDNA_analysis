from types import SimpleNamespace
import gzip
from qc_analysis.scripts.build_trna_position_index import build

def test_existing_output_index_build(tmp_path):
    fasta=tmp_path/'r.fa';fasta.write_text('>chrM\nACGT\n')
    out=tmp_path/'r.out';out.write_text('chrM 1 1 4 Phe GAA 0 0 20.0\n')
    ss=tmp_path/'r.ss';ss.write_text('chrM.trna1 (1-4)\nSeq: ACGU\nStr: ><..\n')
    dest=tmp_path/'r.tsv.gz';a=SimpleNamespace(reference_key='R',fasta=str(fasta),trnascan_out=str(out),trnascan_ss=str(ss),run_trnascan=False,output=str(dest),overwrite=False,chrom_normalization='none',max_sequence_mismatch_rate=0)
    row=build(a);assert row['n_index_rows']==4 and row['n_positive_strand_trna']==1
    with gzip.open(dest,'rt') as h:assert 'paired_base_rna' in h.readline()
