import subprocess, sys
from pathlib import Path

ROOT=Path(__file__).parents[2]
def test_multiallelic_builder(tmp_path):
    v=tmp_path/'S1.vcf'; v.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchrM\t4\t.\tA\tG,T\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:80,10,10:100:0.1,0.1\n')
    m=tmp_path/'metadata.tsv';m.write_text('sample\tspecies\nS1\tSpecies_one\n')
    o=tmp_path/'out.tsv'
    subprocess.run([sys.executable,str(ROOT/'qc_analysis/scripts/build_intraspecies_variant_table.py'),'--vcf-dir',str(tmp_path),'--metadata',str(m),'--output',str(o),'--pass-only','--min-dp','100'],check=True)
    lines=o.read_text().splitlines(); assert len(lines)==3 and lines[1].split('\t')[-1]=='10'
