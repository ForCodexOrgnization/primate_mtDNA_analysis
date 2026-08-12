import csv, gzip, subprocess, sys, shutil, pytest
from pathlib import Path

ROOT=Path(__file__).parents[2]
INTRA=ROOT/"qc_analysis/scripts/run_intraspecies_contamination.py"
FINAL=ROOT/"qc_analysis/scripts/run_final_filter.py"
sys.path.insert(0,str(ROOT))
from qc_analysis.scripts.run_final_filter import find_vcf

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
    if not ((shutil.which("bgzip") and shutil.which("tabix")) or __import__("importlib").util.find_spec("pysam")):
        pytest.skip("pysam or bgzip/tabix is required for production VCF output")
    collected=tmp_path/"collected";(collected/"collected_vcf").mkdir(parents=True)
    (collected/"collected_cov").mkdir();(collected/"collected_mtcn").mkdir();(collected/"reports").mkdir()
    write_vcf(collected/"collected_vcf/A.vcf","A",[(1,.1),(2,.2)]);write_vcf(collected/"collected_vcf/B.vcf","B",[(1,.1)])
    (collected/"collected_cov/A.cov.tsv").write_text("x\n");(collected/"collected_mtcn/A.mtcn.tsv").write_text("x\n")
    (collected/"reports/variant_calling_collection_summary.tsv").write_text("sample\tspecies\nA\tsp\nB\tsp\n")
    intra=tmp_path/"intra.tsv";intra.write_text("sample\tcontamination_status\nA\tno_strong_evidence\nB\thigh_confidence_contaminated\n")
    sample_qc=tmp_path/"sample_qc.tsv";sample_qc.write_text("sample\tqc_status\tfailed_criteria\nA\tPASS\t\nB\tPASS\t\n")
    downstream=tmp_path/"rrna";downstream.mkdir();write_vcf(downstream/"A.lifted.rrna.vcf","A",[(1,.1),(2,.2)]);write_vcf(downstream/"B.lifted.rrna.vcf","B",[(1,.1)])
    flags=tmp_path/"flags.tsv";flags.write_text("sample\tCHROM\tPOS\tREF\tALT\tqc_status\nA\tchrM\t2\tA\tG\tFAIL\n")
    out=tmp_path/"final";cfg=tmp_path/"qc.yaml";cfg.write_text(f"""final_filter:
  enabled: true
  collected_dir: {collected}
  output_dir: {out}
    sample_reports:
    intraspecies:
      path: {intra}
    sample_qc:
      path: {sample_qc}
  vcf_sources: {downstream}
  variant_reports:
    variant_qc:
      path: {flags}
      coordinate_system: human
""")
    run=subprocess.run([sys.executable,str(FINAL),"--config",str(cfg)],text=True,capture_output=True)
    assert run.returncode==0,run.stderr
    assert (out/"final_vcf/A.final.vcf.gz").is_file() and (out/"final_vcf/A.final.vcf.gz.tbi").is_file()
    with gzip.open(out/"final_vcf/A.final.vcf.gz","rt") as h: assert "\t2\t" not in h.read()
    assert (out/"final_cov/A.cov.tsv").is_file() and (out/"final_mtcn/A.mtcn.tsv").is_file()
    for name in ("final_sample_qc.tsv","final_variant_qc.tsv","final_filter_summary.tsv"):assert (out/"reports"/name).is_file()
    with (out/"reports/final_variant_qc.tsv").open() as h: variants=list(csv.DictReader(h,delimiter="\t"))
    assert variants[0]["human_chrom"]=="chrM"
    assert variants[0]["source_chrom"]==variants[0]["original_chrom"]=="NOT_AVAILABLE"

def test_exact_vcf_resolution_does_not_confuse_sample_prefixes(tmp_path):
    write_vcf(tmp_path/"ABC10.lifted.rrna.vcf","ABC10",[(1,.1)])
    assert find_vcf(tmp_path,"ABC1","{sample}.lifted.rrna.vcf") is None
    write_vcf(tmp_path/"ABC1.lifted.rrna.vcf","ABC1",[(1,.1)])
    assert find_vcf(tmp_path,"ABC1","{sample}.lifted.rrna.vcf").name=="ABC1.lifted.rrna.vcf"

def test_exact_vcf_resolution_rejects_compressed_and_plain_ambiguity(tmp_path):
    path=tmp_path/"ABC1.lifted.rrna.vcf";write_vcf(path,"ABC1",[(1,.1)])
    with gzip.open(str(path)+".gz","wt") as h:h.write(path.read_text())
    with pytest.raises(ValueError,match="ambiguous VCF source"):find_vcf(tmp_path,"ABC1","{sample}.lifted.rrna.vcf")

def test_original_coordinate_variant_report_is_rejected(tmp_path):
    # A present report with generic coordinates must explicitly declare that
    # those coordinates are post-liftover human coordinates.
    collected=tmp_path/'collected';(collected/'reports').mkdir(parents=True)
    (collected/'reports/variant_calling_collection_summary.tsv').write_text('sample\tspecies\nA\tsp\n')
    intra=tmp_path/'intra.tsv';intra.write_text('sample\tcontamination_status\nA\tPASS\n')
    sq=tmp_path/'sq.tsv';sq.write_text('sample\tqc_status\nA\tPASS\n')
    flags=tmp_path/'flags.tsv';flags.write_text('sample\tCHROM\tPOS\tREF\tALT\tstatus\nA\tMT\t1\tA\tG\tFAIL\n')
    cfg=tmp_path/'c.yaml';cfg.write_text(f'''final_filter:
  collected_dir: {collected}
  output_dir: {tmp_path}/out
  sample_reports:
    intraspecies:
      path: {intra}
    sample_qc:
      path: {sq}
  variant_reports:
    original_qc:
      path: {flags}
      coordinate_system: source
''')
    run=subprocess.run([sys.executable,str(FINAL),'--config',str(cfg)],text=True,capture_output=True)
    assert run.returncode==2 and 'coordinate_system is unknown or incompatible' in run.stderr

def test_final_filter_fails_when_required_report_missing(tmp_path):
    collected=tmp_path/'collected';(collected/'reports').mkdir(parents=True)
    (collected/'reports/variant_calling_collection_summary.tsv').write_text('sample\tspecies\nA\tsp\n')
    cfg=tmp_path/'c.yaml';cfg.write_text(f'final_filter:\n  collected_dir: {collected}\n  output_dir: {tmp_path}/out\n')
    run=subprocess.run([sys.executable,str(FINAL),'--config',str(cfg)],text=True,capture_output=True)
    assert run.returncode!=0 and 'missing required sample report' in run.stderr
