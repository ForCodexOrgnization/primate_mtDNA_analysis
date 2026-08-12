import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).parents[2];SCRIPT=ROOT/'qc_analysis/scripts/run_sample_variant_filtering.py'
def test_five_qc_boundaries(tmp_path):
 inp=tmp_path/'summary.tsv';inp.write_text('sample\tspecies\tmt_median_coverage\tPercent_100\tnuclear_median_coverage\tmtcn_median\tMAD\nPASS\tsp\t100\t90\t20\t40\t0.4999\nMAD_BOUNDARY\tsp\t100\t90\t20\t40\t0.5\n')
 out=tmp_path/'out';cfg=tmp_path/'c.yaml';cfg.write_text(f'sample_variant_filtering:\n  enabled: true\n  input_summary: {inp}\n  output_dir: {out}\n')
 p=subprocess.run([sys.executable,str(SCRIPT),'--config',str(cfg)],capture_output=True,text=True);assert p.returncode==0,p.stderr
 with (out/'reports/sample_qc.tsv').open() as h:rows={r['sample']:r for r in csv.DictReader(h,delimiter='\t')}
 assert rows['PASS']['qc_status']=='PASS';assert rows['MAD_BOUNDARY']['qc_status']=='FAIL';assert rows['MAD_BOUNDARY']['failed_criteria']=='high_MAD'
