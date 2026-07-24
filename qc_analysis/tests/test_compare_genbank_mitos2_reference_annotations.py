import csv, tempfile, unittest
from pathlib import Path
from qc_analysis.scripts.compare_genbank_mitos2_reference_annotations import compare
from qc_analysis.lib.reference_utils import normalized_fasta_sequence_sha256

class CompareReferenceAnnotationsTests(unittest.TestCase):
 def rows(self,key,sha,gene='MT-ND1',end=3,strand='+'):
  return [dict(reference_key=key,coordinate_reference_sequence_sha256=sha,coordinate_reference_accession='A.1',gene=gene,pos=str(i),ref_base_genome='ATG'[i-1],strand=strand,codon_index='1',codon_pos_in_triplet=str(i),codon_seq='ATG') for i in range(1,end+1)]
 def write(self,p,rows):
  with p.open('w',newline='') as h:
   w=csv.DictWriter(h,fieldnames=sorted({k for r in rows for k in r}),delimiter='\t');w.writeheader();w.writerows(rows)
 def test_exact_and_sequence_mismatch_are_separated(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td); fa=d/'a.fa';fa.write_text('>x\nATG\n'); sha=normalized_fasta_sequence_sha256(fa)['sequence_sha256']
   g,m=d/'g.tsv',d/'m.tsv'; self.write(g,self.rows('g',sha));self.write(m,self.rows('m',sha))
   out,summ,diag=d/'o.tsv',d/'s.tsv',d/'d.tsv';compare(g,m,out,summ,diag)
   with out.open() as h: self.assertEqual(next(csv.DictReader(h,delimiter='\t'))['gene_comparison_category'],'exact_match')
   with summ.open() as h: self.assertEqual(next(csv.DictReader(h,delimiter='\t'))['reference_comparison_category'],'gene_missing') # only one of 13 present
   self.write(m,self.rows('m','f'*64)); compare(g,m,out,summ,diag,fail_no_shared=False)
   with diag.open() as h: self.assertEqual(next(csv.DictReader(h,delimiter='\t'))['sequence_compatibility_category'],'sequence_mismatch')
   with out.open() as h: self.assertEqual(list(csv.DictReader(h,delimiter='\t')),[])
