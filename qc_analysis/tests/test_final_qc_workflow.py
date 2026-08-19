import csv, gzip, subprocess, sys, shutil, pytest
from pathlib import Path

ROOT=Path(__file__).parents[2]
INTRA=ROOT/"qc_analysis/scripts/run_intraspecies_contamination.py"
FINAL=ROOT/"qc_analysis/scripts/run_final_filter.py"
sys.path.insert(0,str(ROOT))
from qc_analysis.scripts.run_final_filter import call_class, find_vcf, sort_plain_vcf, variant_classes

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
    downstream=tmp_path/"rrna";downstream.mkdir();write_vcf(downstream/"A.lifted.codon.trna.rrna.vcf","A",[(1,.1),(2,.2)]);write_vcf(downstream/"B.lifted.codon.trna.rrna.vcf","B",[(1,.1)])
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

def final_filter_fixture(tmp_path, records, filtered_positions=()):
    collected=tmp_path/"collected";(collected/"reports").mkdir(parents=True)
    (collected/"reports/variant_calling_collection_summary.tsv").write_text("sample\tspecies\nA\tsp\n")
    intra=tmp_path/"intra.tsv";intra.write_text("sample\tcontamination_status\nA\tno_strong_evidence\n")
    sample_qc=tmp_path/"sample_qc.tsv";sample_qc.write_text("sample\tqc_status\nA\tPASS\n")
    downstream=tmp_path/"rrna";downstream.mkdir();write_vcf(downstream/"A.lifted.codon.trna.rrna.vcf","A",records)
    flags=tmp_path/"flags.tsv"
    flags.write_text(
        "sample\tCHROM\tPOS\tREF\tALT\tqc_status\n"
        + "".join(f"A\tchrM\t{pos}\tA\tG\tFAIL\n" for pos in filtered_positions)
    )
    out=tmp_path/"final";cfg=tmp_path/"qc.yaml"
    cfg.write_text(f"""final_filter:
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
    return out

def final_vcf_positions(path):
    with gzip.open(path,"rt") as handle:
        return [int(line.split("\t")[1]) for line in handle if not line.startswith("#")]

def test_final_vcf_is_coordinate_sorted_before_indexing(tmp_path):
    if not ((shutil.which("bgzip") and shutil.which("tabix")) or __import__("importlib").util.find_spec("pysam")):
        pytest.skip("pysam or bgzip/tabix is required for production VCF output")
    out=final_filter_fixture(tmp_path,[(16465,.1),(74,.2),(150,.3)])
    final=out/"final_vcf/A.final.vcf.gz"
    assert final_vcf_positions(final)==[74,150,16465]
    assert final.is_file()
    assert Path(str(final)+".tbi").is_file()

def test_final_vcf_sorting_occurs_after_variant_filtering(tmp_path):
    if not ((shutil.which("bgzip") and shutil.which("tabix")) or __import__("importlib").util.find_spec("pysam")):
        pytest.skip("pysam or bgzip/tabix is required for production VCF output")
    out=final_filter_fixture(tmp_path,[(16465,.1),(74,.2),(150,.3)],filtered_positions=(150,))
    assert final_vcf_positions(out/"final_vcf/A.final.vcf.gz")==[74,16465]
    with (out/"reports/final_variant_qc.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert [row["human_pos"] for row in rows]==["16465","74","150"]
    assert [row["final_variant_status"] for row in rows]==["PASS","PASS","FAIL"]

def test_sort_plain_vcf_sorts_positions_numerically(tmp_path):
    source=tmp_path/"input.vcf";output=tmp_path/"sorted.vcf"
    write_vcf(source,"A",[(1000,.1),(100,.1),(9,.1),(74,.1)])
    original_headers=[line for line in source.read_text().splitlines(keepends=True) if line.startswith("#")]
    sort_plain_vcf(source,output)
    lines=output.read_text().splitlines(keepends=True)
    assert [line for line in lines if line.startswith("#")]==original_headers
    assert [int(line.split("\t")[1]) for line in lines if not line.startswith("#")]==[9,74,100,1000]

def test_variant_report_classification_boundaries():
    assert [call_class(value) for value in (.95,.10,.099,None)]==["homoplasmic","heteroplasmic","low_af","UNKNOWN"]
    assert variant_classes("A","G")==("SNV","SNV_transition")
    assert variant_classes("C","A")==("SNV","SNV_transversion")
    assert variant_classes("A","AT")==("INDEL","indel")
    assert variant_classes("AC","GT")==("OTHER","other")

def test_final_variant_report_exposes_vcf_sample_and_orthology_annotations(tmp_path):
    if not ((shutil.which("bgzip") and shutil.which("tabix")) or __import__("importlib").util.find_spec("pysam")):
        pytest.skip("pysam or bgzip/tabix is required for production VCF output")
    collected=tmp_path/"collected";(collected/"reports").mkdir(parents=True)
    (collected/"reports/variant_calling_collection_summary.tsv").write_text("sample\tspecies\nA\tPan_troglodytes\n")
    intra=tmp_path/"intra.tsv";intra.write_text("sample\tcontamination_status\nA\tno_strong_evidence\n")
    sample_qc=tmp_path/"sample_qc.tsv"
    sample_qc.write_text("sample\tmt_median_coverage\tPercent_100\tnuclear_median_coverage\tmtcn_median\tMAD\tfailed_criteria\tqc_status\nA\t500\t99.5\t30\t72\t0.1\t\tPASS\n")
    downstream=tmp_path/"rrna";downstream.mkdir();vcf=downstream/"A.lifted.codon.trna.rrna.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
        "chrM\t100\t.\tA\tG\t.\tPASS\tDP=999;MTLIFT_ORIG_CHROM=species;MTLIFT_ORIG_POS=10;MTLIFT_ORIG_REF=C;MTLIFT_ORIG_ALT=T;MTCODON_STATUS=PASS;MTTRNA_STATUS=NO_SPECIES_OR_HUMAN_TRNA;MTRRNA_STATUS=OK\tGT:DP:AF:AD\t0/1:321:0.96:50,50\n"
        "chrM\t101\t.\tC\tA\t.\tq10\tSRC_CHROM=species;SRC_POS=11;SRC_REF=G;SRC_ALT=T\tGT:DP:AF\t0/1:200:0.50\n"
        "chrM\t102\t.\tT\tC\t.\tPASS\t.\tGT:DP:AF\t0/1:150:0.05\n"
        "chrM\t103\t.\tA\tAT\t.\tPASS\tDP=44\tGT:AD\t0/1:10,90\n"
    )
    orthology=tmp_path/"orthology.tsv"
    orthology.write_text(
        "sample\thuman_chrom\thuman_pos\thuman_ref\thuman_alt\tregion_type\torthology_match_status\torthology_fail_reason\n"
        "A\tchrM\t100\tA\tG\tCDS\tPASS\t\n"
        "A\tchrM\t101\tC\tA\ttRNA\tFAIL\tMISMATCH\n"
    )
    out=tmp_path/"final";cfg=tmp_path/"qc.yaml"
    cfg.write_text(f"""final_filter:
  collected_dir: {collected}
  output_dir: {out}
  sample_reports:
    intraspecies:
      path: {intra}
    sample_qc:
      path: {sample_qc}
  vcf_sources: {downstream}
  variant_reports:
    orthology:
      path: {orthology}
      coordinate_system: human
      status_columns: orthology_match_status
      fail_status: FAIL
""")
    run=subprocess.run([sys.executable,str(FINAL),"--config",str(cfg)],text=True,capture_output=True)
    assert run.returncode==0,run.stderr
    with (out/"reports/final_variant_qc.tsv").open() as handle:
        rows={row["human_pos"]:row for row in csv.DictReader(handle,delimiter="\t")}
    assert rows["100"]["AF"]=="0.96" and rows["100"]["DP"]=="321.0"
    assert rows["103"]["AF"]=="0.9" and rows["103"]["DP"]=="44.0"
    assert [rows[pos]["call_class"] for pos in ("100","101","102","103")]==["homoplasmic","heteroplasmic","low_af","heteroplasmic"]
    assert rows["100"]["snv_type"]=="SNV_transition"
    assert rows["101"]["snv_type"]=="SNV_transversion"
    assert rows["103"]["variant_class"]=="INDEL" and rows["103"]["snv_type"]=="indel"
    assert rows["101"]["vcf_filter"]=="q10"
    assert [rows["100"][field] for field in ("mt_median_coverage","Percent_100","nuclear_median_coverage","mtcn_median","MAD")]==["500","99.5","30","72","0.1"]
    assert rows["100"]["species"]=="Pan_troglodytes" and rows["100"]["intraspecies_status"]=="no_strong_evidence"
    assert rows["100"]["human_contamination_status"]==rows["100"]["interspecies_status"]=="NOT_AVAILABLE"
    assert [rows["100"][field] for field in ("source_chrom","source_pos","source_ref","source_alt")]==["species","10","C","T"]
    assert [rows["101"][field] for field in ("source_chrom","source_pos","source_ref","source_alt")]==["species","11","G","T"]
    assert [rows["100"][field] for field in ("codon_match_status","trna_match_status","rrna_match_status")]==["PASS","NO_SPECIES_OR_HUMAN_TRNA","OK"]
    assert [rows["100"][field] for field in ("region_type","orthology_match_status","orthology_fail_reason")]==["CDS","PASS","NOT_AVAILABLE"]
    assert rows["101"]["final_variant_status"]=="FAIL" and rows["101"]["final_variant_fail_reasons"]=="orthology:FAIL"
    assert final_vcf_positions(out/"final_vcf/A.final.vcf.gz")==[100,102,103]

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
