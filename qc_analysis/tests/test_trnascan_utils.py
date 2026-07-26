import csv,gzip
from pathlib import Path
import pytest
from qc_analysis.lib.trnascan_utils import *
from qc_analysis.lib import trnascan_utils

def inputs(tmp_path,negative=False,sequence='ACGT',structure='><..'):
    begin,end=(4,1) if negative else (1,4)
    out=tmp_path/'x.out';out.write_text(f'Sequence\ttRNA #\tBegin\tEnd\tType\tCodon\tIntron Begin\tIntron End\tScore\nchrM\t1\t{begin}\t{end}\tPhe\tGAA\t0\t0\t42.0\n')
    ss=tmp_path/'x.ss';ss.write_text(f'chrM.trna1 ({begin}-{end})\nType: Phe Anticodon: GAA at 2-4\nSeq: {sequence}\nStr: {structure}\n')
    return out,ss

def test_parse_positive_and_structure(tmp_path):
    out,ss=inputs(tmp_path);r=merge_trnascan_records(parse_trnascan_out(out),parse_trnascan_ss(ss))[0]
    assert r.strand=='+' and r.pairs=={1:2,2:1} and r.genomic_pos(3)==3
    assert infer_structural_elements(r.structure)[1]=='acceptor_stem'

def test_negative_coordinates_and_rna_pairs(tmp_path):
    out,ss=inputs(tmp_path,True,sequence='ACGU')
    fasta=tmp_path/'x.fa';fasta.write_text('>chrM\nACGT\n');dest=tmp_path/'i.tsv.gz'
    result=build_trna_position_index('R',fasta,out,ss,dest,mismatch_rate_threshold=1)
    assert result['rows'][0]['pos']==4 and result['rows'][0]['base_rna']=='A'
    assert result['rows'][0]['paired_genomic_pos']==3
    with gzip.open(dest,'rt') as h:assert next(csv.DictReader(h,delimiter='\t'))['index_format_version']=='2'

@pytest.mark.parametrize('a,b,kind',[('A','U','WC'),('G','U','GU_wobble'),('A','C','non_WC')])
def test_pair_types(a,b,kind):assert trnascan_utils._pair_type(a,b)==kind

def test_malformed_and_no_trna(tmp_path):
    ss=tmp_path/'bad.ss';ss.write_text('nothing\n')
    with pytest.raises(ValueError):parse_trnascan_ss(ss)
    out=tmp_path/'empty.out';out.write_text('Sequence tRNA # Begin End Type Codon X X Score\n')
    good=tmp_path/'good.ss';good.write_text('chrM.trna1 (1-2)\nSeq: AU\nStr: ><\n')
    with pytest.raises(ValueError,match='No tRNAs'):merge_trnascan_records(parse_trnascan_out(out),parse_trnascan_ss(good))

def test_sequence_validation(tmp_path):
    out,ss=inputs(tmp_path,sequence='AAAA');f=tmp_path/'x.fa';f.write_text('>chrM\nACGT\n')
    with pytest.raises(ValueError,match='mismatch rate'):build_trna_position_index('R',f,out,ss,tmp_path/'x.tsv')
