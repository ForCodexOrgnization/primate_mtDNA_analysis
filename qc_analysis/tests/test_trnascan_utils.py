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

def test_iupac_reference_build_preserves_symbols_and_ambiguous_pair(tmp_path):
    out,ss=inputs(tmp_path,sequence='AGAU',structure='><><')
    fasta=tmp_path/'iupac.fa'; fasta.write_text('>chrM\nMRWY\n')
    dest=tmp_path/'iupac.tsv'
    result=build_trna_position_index('R',fasta,out,ss,dest)
    assert result['n_fasta_sequence_mismatch']==0
    assert ''.join(row['base_genomic'] for row in result['rows'])=='MRWY'
    assert result['rows'][0]['pair_status']=='paired'
    assert result['rows'][0]['pair_type']=='ambiguous'
    assert result['rows'][0]['pair_state']=='NA'
    assert validate_trna_index(dest,'R')['n_rows']==4

def test_negative_iupac_reference_orientation(tmp_path):
    out,ss=inputs(tmp_path,negative=True,sequence='WRYK',structure='....')
    fasta=tmp_path/'iupac.fa'; fasta.write_text('>chrM\nMRWY\n')
    result=build_trna_position_index('R',fasta,out,ss,tmp_path/'iupac.tsv')
    assert [row['base_genomic'] for row in result['rows']]==list('YWRM')
    assert [row['base_rna'] for row in result['rows']]==list('RWYK')

def test_malformed_and_no_trna(tmp_path):
    ss=tmp_path/'bad.ss';ss.write_text('nothing\n')
    with pytest.raises(ValueError):parse_trnascan_ss(ss)
    out=tmp_path/'empty.out';out.write_text('Sequence tRNA # Begin End Type Codon X X Score\n')
    good=tmp_path/'good.ss';good.write_text('chrM.trna1 (1-2)\nSeq: AU\nStr: ><\n')
    with pytest.raises(ValueError,match='No tRNAs'):merge_trnascan_records(parse_trnascan_out(out),parse_trnascan_ss(good))

def test_sequence_validation(tmp_path):
    out,ss=inputs(tmp_path,sequence='AAAA');f=tmp_path/'x.fa';f.write_text('>chrM\nACGT\n')
    with pytest.raises(ValueError,match='mismatch rate'):build_trna_position_index('R',f,out,ss,tmp_path/'x.tsv')

def test_reordered_ss_records_match_by_id(tmp_path):
    out=tmp_path/'x.out'; out.write_text('chrM 1 1 4 Phe GAA 0 0 20\nchrM 2 5 8 Leu TAA 0 0 20\n')
    ss=tmp_path/'x.ss'; ss.write_text('chrM.trna2 (5-8)\nType: Leu Anticodon: TAA\nSeq: AAAA\nStr: ....\nchrM.trna1 (1-4)\nType: Phe Anticodon: GAA\nSeq: CCCC\nStr: ....\n')
    records=merge_trnascan_records(parse_trnascan_out(out),parse_trnascan_ss(ss))
    assert [r.sequence for r in records]==['CCCC','AAAA']

def test_topology_structural_inference_and_noncanonical_arm():
    canonical='>>>>...<<<<'  # topology is used even though canonical coordinates are absent
    labels=infer_structural_elements(canonical)
    assert labels[1]=='acceptor_stem' and labels[5]=='connector'
    shortened='>>..<<'
    assert infer_structural_elements(shortened)[1]=='acceptor_stem'
    assert all(infer_structural_elements('......').values())  # documented canonical fallback

def _successful_trnascan(run_args, **kwargs):
    for option in ('-o', '-f'):
        Path(run_args[run_args.index(option) + 1]).write_text('new output\n')

def test_run_trnascan_removes_exact_prefix_outputs_and_is_noninteractive(tmp_path, monkeypatch):
    prefix=tmp_path/'reference'; other=tmp_path/'reference_extra.trnascan.out'
    stale=[Path(str(prefix)+suffix) for suffix in
           ('.trnascan.out','.trnascan.ss','.trnascan.stats','.trnascan.bed','.trnascan.fa')]
    for path in stale:path.write_text('stale')
    other.write_text('unrelated')
    observed={}
    def fake_run(command, **kwargs):
        observed.update(command=command,kwargs=kwargs)
        assert all(not path.exists() for path in stale)
        _successful_trnascan(command, **kwargs)
    monkeypatch.setattr(subprocess,'run',fake_run)

    made=run_trnascan(tmp_path/'input.fa',prefix)

    assert other.read_text()=='unrelated'
    assert '--forceow' in observed['command']
    assert 'input' not in observed['kwargs'] and 'stdin' not in observed['kwargs']
    assert observed['kwargs']['check'] is True and observed['kwargs']['timeout']==3600
    assert observed['kwargs']['stdout'].name==str(prefix)+'.stdout.log'
    assert observed['kwargs']['stderr'].name==str(prefix)+'.stderr.log'
    assert made['out'].read_text()=='new output\n'

def test_run_trnascan_clean_directory_is_unchanged(tmp_path, monkeypatch):
    prefix=tmp_path/'clean'
    monkeypatch.setattr(subprocess,'run',_successful_trnascan)
    result=run_trnascan(tmp_path/'input.fa',prefix)
    assert result['out'].is_file() and result['ss'].is_file()

def test_run_trnascan_without_overwrite_refuses_to_launch(tmp_path, monkeypatch):
    prefix=tmp_path/'reference'; Path(str(prefix)+'.trnascan.out').write_text('stale')
    called=False
    def fake_run(*args,**kwargs):
        nonlocal called;called=True
    monkeypatch.setattr(subprocess,'run',fake_run)
    with pytest.raises(FileExistsError,match='enable overwrite'):
        run_trnascan(tmp_path/'input.fa',prefix,overwrite=False)
    assert not called

def test_run_trnascan_failure_is_captured_in_per_reference_log(tmp_path, monkeypatch):
    prefix=tmp_path/'failed'
    def fail(command,stdout,stderr,**kwargs):
        stderr.write('diagnostic from tool\n');stderr.flush()
        raise subprocess.CalledProcessError(7,command)
    monkeypatch.setattr(subprocess,'run',fail)
    with pytest.raises(RuntimeError,match=r'exit code 7.*failed\.stderr\.log'):
        run_trnascan(tmp_path/'input.fa',prefix)
    assert Path(str(prefix)+'.stderr.log').read_text()=='diagnostic from tool\n'
