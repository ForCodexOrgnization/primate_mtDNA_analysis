import csv, gzip, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).parents[2]
INTRA=ROOT/"qc_analysis/scripts/run_intraspecies_contamination.py"
FINAL=ROOT/"qc_analysis/scripts/run_final_filter.py"

def write_vcf(path, sample, records):
    text="##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample+"\n"
    for pos,af in records:text+=f"chrM\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT:DP:AF\t0/1:200:{af}\n"
    path.write_text(text)

def test_python_intraspecies_report_contract(tmp_path):
    vcf=tmp_path/"vcf";vcf.mkdir();write_vcf(vcf/"A.vcf","A",[(1,.1),(2,.1),(3,.1),(4,.1),(5,.1)])
    write_vcf(vcf/"B.vcf","B",[(1,.999),(2,.999),(3,.999),(4,.999),(5,.999)])
    summary=tmp_path/"summary.tsv";summary.write_text("sample\tspecies\nA\tsp\nB\tsp\n")
    out=tmp_path/"out";cfg=tmp_path/"qc.yaml"
    cfg.write_text(f"""intraspecies_contamination:
  enabled: true
  vcf_dir: {vcf}
  sample_summary: {summary}
  outdir: {out}
  overwrite: true
  dp_min: 100
  use_snv_only: true
  pass_only: true
""")
    run=subprocess.run([sys.executable,str(INTRA),"--config",str(cfg)],text=True,capture_output=True)
    assert run.returncode==0,run.stderr
    report=out/"reports/intraspecies_contamination_report.tsv"
    with report.open() as h: rows=list(csv.DictReader(h,delimiter="\t"))
    assert {r["sample"] for r in rows}=={"A","B"}
    assert {"contamination_status","qc_status","n_depressed_anchor"} <= rows[0].keys()
    assert (out/"logs/intraspecies_contamination.log").is_file()

def test_final_filter_excludes_failed_sample_and_failed_variant(tmp_path):
    collected=tmp_path/"collected";(collected/"collected_vcf").mkdir(parents=True)
    (collected/"collected_cov").mkdir();(collected/"collected_mtcn").mkdir();(collected/"reports").mkdir()
    write_vcf(collected/"collected_vcf/A.vcf","A",[(1,.1),(2,.2)]);write_vcf(collected/"collected_vcf/B.vcf","B",[(1,.1)])
    (collected/"collected_cov/A.cov.tsv").write_text("x\n");(collected/"collected_mtcn/A.mtcn.tsv").write_text("x\n")
    (collected/"reports/variant_calling_collection_summary.tsv").write_text("sample\tspecies\nA\tsp\nB\tsp\n")
    intra=tmp_path/"intra.tsv";intra.write_text("sample\tcontamination_status\nA\tno_strong_evidence\nB\thigh_confidence_contaminated\n")
    flags=tmp_path/"flags.tsv";flags.write_text("sample\tCHROM\tPOS\tREF\tALT\tqc_status\nA\tchrM\t2\tA\tG\tFAIL\n")
    out=tmp_path/"final";cfg=tmp_path/"qc.yaml";cfg.write_text(f"""final_filter:
  enabled: true
  collected_dir: {collected}
  output_dir: {out}
  sample_reports:
    intraspecies:
      path: {intra}
  variant_reports:
    variant_qc:
      path: {flags}
""")
    run=subprocess.run([sys.executable,str(FINAL),"--config",str(cfg)],text=True,capture_output=True)
    assert run.returncode==0,run.stderr
    assert (out/"final_vcf/A.vcf").is_file() and not (out/"final_vcf/B.vcf").exists()
    assert "\t2\t" not in (out/"final_vcf/A.vcf").read_text()
    assert (out/"final_cov/A.cov.tsv").is_file() and (out/"final_mtcn/A.mtcn.tsv").is_file()
    for name in ("final_sample_qc.tsv","final_variant_qc.tsv","final_filter_summary.tsv"):assert (out/"reports"/name).is_file()
