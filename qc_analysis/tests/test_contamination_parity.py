import csv,math,shutil,subprocess,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).parents[2];VAR=ROOT/'qc_analysis/tests/fixtures/contamination_parity_variants.tsv';GOLD=ROOT/'qc_analysis/tests/fixtures/contamination_parity_golden.tsv'
FIELDS='sample species n_lowA best_source_sample best_overlap best_frac_lowA_in_highB mt_high_hets_contamination mt_high_hets_mode n_mirror_pairs n_low_variants_with_mirror mirror_low_fraction contamination_flag_candidate contamination_flag_highconf contamination_status'.split()
def rows(p):
 with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
def same(a,b):
 for f in FIELDS:
  av,bv=a.get(f,a.get(f.capitalize(),'')),b.get(f,b.get(f.capitalize(),''))
  try:
   if math.isnan(float(av)) and math.isnan(float(bv)):continue
   assert float(av)==pytest.approx(float(bv),abs=1e-10),(f,av,bv)
  except ValueError:assert str(av).lower()==str(bv).lower(),(f,av,bv)
def test_python_matches_r_generated_golden(tmp_path):
 out=tmp_path/'py';p=subprocess.run([sys.executable,str(ROOT/'qc_analysis/scripts/run_intraspecies_contamination.py'),'--variant-table',str(VAR),'--outdir',str(out),'--overwrite'],capture_output=True,text=True);assert p.returncode==0,p.stderr
 actual=rows(out/'reports/intraspecies_contamination_report.tsv');gold=rows(GOLD);assert len(actual)==len(gold)
 for a,b in zip(actual,gold):same(a,b)
def test_r_and_python_reference_parity(tmp_path):
 if not shutil.which('Rscript'):pytest.skip('Rscript unavailable; Python golden parity remains enforced')
 py=tmp_path/'py';rr=tmp_path/'r';subprocess.run([sys.executable,str(ROOT/'qc_analysis/scripts/run_intraspecies_contamination.py'),'--variant-table',str(VAR),'--outdir',str(py),'--overwrite'],check=True)
 subprocess.run(['Rscript',str(ROOT/'qc_analysis/validation/contamination_reference.R'),'--variant-table',str(VAR),'--outdir',str(rr),'--overwrite'],check=True,cwd=ROOT)
 pr=rows(py/'reports/intraspecies_contamination_report.tsv');rraw=rows(rr/'tables/final_contamination_summary.tsv');mapped=[]
 for x in rraw:
  y={k.lower():v for k,v in x.items()};y['sample']=x['Sample'];y['species']=x['Species'];mapped.append(y)
 for a,b in zip(sorted(pr,key=lambda x:x['sample']),sorted(mapped,key=lambda x:x['sample'])):same(a,b)
