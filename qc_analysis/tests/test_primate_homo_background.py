import csv, importlib.util, json
from pathlib import Path
import pytest

SCRIPT=Path(__file__).parents[1]/"scripts/build_primate_homo_background.py"
spec=importlib.util.spec_from_file_location("background",SCRIPT);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def test_metadata_reads_headered_manifest_case_insensitively(tmp_path):
    path=tmp_path/"samples.tsv"
    path.write_text("# cohort metadata\n\nSAMPLE\tSpecies\tGenus\tFamily\tnote\nA\tPan_troglodytes\tPan\tHominidae\talpha\nB\tGorilla_gorilla\tGorilla\tHominidae\tbeta\n")
    rows=mod.metadata(path)
    assert set(rows)=={"A","B"}
    assert rows["A"]=={"sample":"A","species":"Pan_troglodytes","genus":"Pan","family":"Hominidae","note":"alpha"}

def test_metadata_reads_legacy_headerless_manifest(tmp_path):
    path=tmp_path/"samples.tsv"
    path.write_text("# sample and species\n\nA\tPan_troglodytes\textra\nB\tGorilla_gorilla\n")
    rows=mod.metadata(path)
    assert set(rows)=={"A","B"}
    assert rows["A"]=={"sample":"A","species":"Pan_troglodytes"}
    assert rows["B"]=={"sample":"B","species":"Gorilla_gorilla"}

def test_metadata_rejects_duplicate_sample_ids(tmp_path):
    path=tmp_path/"samples.tsv"
    path.write_text("sample\tspecies\nA\tPan_troglodytes\nA\tGorilla_gorilla\n")
    with pytest.raises(ValueError,match="Duplicate sample ID 'A'.*lines 2 and 3"):
        mod.metadata(path)

def test_existing_annotation_statuses_are_consolidated_without_new_matching_rules():
    tiers={"HIGH_CONF_STEM","HIGH_CONF_LOOP"}
    assert mod.status({"MTCODON_STATUS":"PASS"},tiers)[:2]==("CDS","PASS")
    assert mod.status({"MTCODON_STATUS":"AMBIGUOUS_CODON"},tiers)[:2]==("CDS","AMBIGUOUS")
    assert mod.status({"MTCODON_STATUS":"SKIPPED_NONCODING","MTTRNA_STATUS":"OK","MTTRNA_STRICT_MATCH":"no"},tiers)[:2]==("tRNA","FAIL")
    assert mod.status({"MTCODON_STATUS":"SKIPPED_NONCODING","MTTRNA_STATUS":"NO_SPECIES_OR_HUMAN_TRNA","MTRRNA_STATUS":"NO_SPECIES_OR_HUMAN_RRNA"},tiers)[:2]==("noncoding","NOT_APPLICABLE")

@pytest.mark.parametrize(("info","expected"),[
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"HIGH_CONF_STEM"},("rRNA","PASS","")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"HIGH_CONF_LOOP"},("rRNA","PASS","")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"LOW_CONF"},("rRNA","FAIL","RRNA_LOW_CONF")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"STRUCTURE_DISCORDANT"},("rRNA","FAIL","RRNA_STRUCTURE_DISCORDANT")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"STRUCTURE_UNKNOWN"},("rRNA","FAIL","RRNA_STRUCTURE_UNKNOWN")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"no","MTRRNA_MATCH_TIER":"HIGH_CONF_STEM"},("rRNA","FAIL","REGION_MATCH_NO")),
    ({"MTRRNA_STATUS":"NO_SPECIES_RRNA"},("rRNA","FAIL","NO_SPECIES_RRNA")),
    ({"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes"},("rRNA","FAIL","RRNA_MATCH_TIER_MISSING")),
])
def test_rrna_status_requires_an_accepted_structural_tier(info,expected):
    assert mod.status(info,{"HIGH_CONF_STEM","HIGH_CONF_LOOP"})==expected

def test_rrna_status_uses_configured_tier_set():
    info={"MTRRNA_STATUS":"OK","MTRRNA_REGION_MATCH":"yes","MTRRNA_MATCH_TIER":"LOW_CONF"}
    assert mod.status(info,{"LOW_CONF"})==("rRNA","PASS","")

def test_empty_accepted_rrna_tier_configuration_fails_clearly(tmp_path,monkeypatch):
    cfg=tmp_path/"c.yaml";cfg.write_text('primate_homo_background:\n  settings:\n    accepted_rrna_match_tiers: ""\n')
    monkeypatch.setattr("sys.argv",[str(SCRIPT),"--config",str(cfg)])
    with pytest.raises(ValueError,match="accepted_rrna_match_tiers must contain at least one tier"):
        mod.main()

def test_background_requires_homoplasmy_depth_filter_snv_and_match(tmp_path,monkeypatch):
    root=tmp_path;vcfs=root/"rrna";vcfs.mkdir();out=root/"out";ortho=root/"ortho"
    meta=root/"samples.tsv";meta.write_text("sample\tspecies\tgenus\tfamily\ns\tsp\tg\tf\n")
    marker=root/"markers.tsv";marker.write_text("human_pos\thuman_ref\thuman_alt\n100\tA\tG\n")
    records=[
      (100,"A","G","PASS","MTCODON_STATUS=PASS",.96,120),
      (101,"A","G","PASS","MTCODON_STATUS=PASS",.94,120),
      (102,"A","G","PASS","MTCODON_STATUS=MISMATCH",.99,120),
      (103,"A","G","q10","MTCODON_STATUS=PASS",.99,120),
      (104,"A","AT","PASS","MTCODON_STATUS=PASS",.99,120),
    ]
    with (vcfs/"s.lifted.codon.trna.rrna.vcf").open("w") as h:
      h.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n")
      for pos,ref,alt,flt,info,af,dp in records:h.write(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t{info}\tAF:DP\t{af}:{dp}\n")
    cfg=root/"c.yaml";cfg.write_text(f'''primate_homo_background:
  paths:
    input_vcf_dir: {vcfs}
    sample_ref_file: {meta}
    output_dir: {out}
    orthology_reports_dir: {ortho}
    human_marker_table: {marker}
  settings:
    homoplasmy_af_min: 0.95
    dp_min: 100
    pass_only: true
    snv_only: true
    accepted_orthology_statuses: PASS
''')
    monkeypatch.setattr("sys.argv",[str(SCRIPT),"--config",str(cfg)]);mod.main()
    with (out/"primate_homo_background.tsv").open() as h:rows=list(csv.DictReader(h,delimiter="\t"))
    assert [(r["human_pos"],r["human_ref"],r["human_alt"]) for r in rows]==[("100","A","G")]
    with (ortho/"orthology_match_report.tsv").open() as h:orthology=list(csv.DictReader(h,delimiter="\t"))
    assert len(orthology)==len(records)  # report retains every input record
    provenance=json.loads((out/"primate_homo_background_metadata.json").read_text())
    assert provenance["accepted_rrna_match_tiers"]==["HIGH_CONF_LOOP","HIGH_CONF_STEM"]

def test_workflow_places_background_before_human():
    text=(Path(__file__).parents[1]/"scripts/run_qc_preprocessing.sh").read_text()
    graph=text[text.index("local steps=("):text.index("for s in",text.index("local steps=("))]
    assert graph.index("rrna_match") < graph.index("build_primate_homo_background") < graph.index("human_contamination") < graph.index("final_filter")
