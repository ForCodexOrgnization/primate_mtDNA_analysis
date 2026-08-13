import importlib.util
from pathlib import Path

SCRIPT=Path(__file__).parents[1]/"scripts/run_human_contamination.py"
spec=importlib.util.spec_from_file_location("human_qc",SCRIPT);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def config(**changes):
    cfg={"variant_filters":{"dp_min":100,"pass_only":True,"snv_only":True,"low_vaf_min":.01,"low_vaf_max":.5},
         "marker_screen":{"min_low_variants_for_screen":6,"min_human_marker_hits":6,"min_fraction_low_variants_human_marker":.6,"include_back_mutations_in_candidate_screen":True,"include_back_mutations_in_fail_screen":False},
         "vaf_coherence":{"enabled":True,"tolerance":.03,"min_fraction_markers_coherent":.7},
         "control_region":{"start":16000,"end":576,"min_non_control_region_hits_for_fail":3},
         "classification":{"require_baseline_marker_screen_for_fail":True,"require_vaf_coherence_for_fail":True,"require_non_control_markers_for_fail":True},
         "haplogrep":{"input_vaf_min":.01,"input_vaf_max":.5,"exclude_back_mutations":False,"require_phylotree_marker":True}}
    for section,value in changes.items():cfg[section].update(value)
    return cfg

def write_vcf(tmp_path, records):
    p=tmp_path/"s.lifted.raw.vcf"
    with p.open("w") as h:
        h.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for pos,alt,af in records:
            value=".:95,5:100" if af is None else f"{af}:95,5:100"
            h.write(f"chrM\t{pos}\t.\tA\t{alt}\t.\tPASS\tSRC_ALT=T;LIFTOVER_ALLELE_STATUS=ALT_REF_FLIP\tAF:AD:DP\t{value}\n")
    return p

def markers(positions, back=()):
    return {(p,"G"):{"marker":f"{p}G"+("!" if p in back else ""),"is_back_mutation":p in back} for p in positions}

def test_statuses_and_final_lifted_allele(tmp_path):
    positions=[700,800,900,1000,1100,1200]
    audit=[];row,_=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(p,"G",x) for p,x in zip(positions,[.047,.051,.050,.044,.053,.048])]),markers(positions),config(),audit)
    assert row["human_contamination_status"]=="FAIL"
    assert row["n_human_marker_hits"]==6 and all(x["liftover_allele_status"]=="ALT_REF_FLIP" for x in audit)
    # SRC_ALT=T cannot prevent matching the canonical lifted ALT=G.

def test_tiny_denominator_and_no_markers(tmp_path):
    row,_=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(700,"G",.05)]),markers([700]),config(),[])
    assert row["human_contamination_status"]=="INSUFFICIENT_DATA"
    row,_=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(x,"T",.05) for x in range(700,1300,100)]),markers([]),config(),[])
    assert row["human_contamination_status"]=="PASS"

def test_incoherent_and_control_only_are_candidates(tmp_path):
    positions=[700,800,900,1000,1100,1200]
    row,_=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(p,"G",x) for p,x in zip(positions,[.01,.1,.2,.3,.4,.5])]),markers(positions),config(),[])
    assert row["human_contamination_status"]=="CANDIDATE"
    positions=[73,146,152,263,300,500]
    row,_=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(p,"G",.05) for p in positions]),markers(positions),config(),[])
    assert row["human_contamination_status"]=="CANDIDATE" and not row["non_control_marker_pass"]

def test_ad_fallback_missing_and_marker_dedup(tmp_path):
    p=write_vcf(tmp_path,[(700,"G",None)])
    variant=list(mod.parse_vcf(p))[0]
    assert variant["af"]==.05 and variant["af_source"]=="CALCULATED_FROM_AD"
    table=tmp_path/"markers.tsv";table.write_text("marker\tpos\talt\tis_back_mutation\n700G\t700\tG\tfalse\n700G\t700\tG\tfalse\nbad-del\t0\tA\tfalse\n")
    loaded,qc=mod.load_markers(table);assert len(loaded)==1 and qc["duplicates"]==1 and qc["excluded"]==1

def test_back_mutation_exclusion_and_marker_only_profile(tmp_path):
    positions=[700,800,900,1000,1100,1200]
    cfg=config(marker_screen={"include_back_mutations_in_candidate_screen":False})
    row,selected=mod.analyze_sample("s","sp",write_vcf(tmp_path,[(p,"G",.05) for p in positions]+[(1300,"T",.05)]),markers(positions,back=[700]),cfg,[])
    assert row["n_human_marker_hits"]==5
    assert all(k in markers(positions,back=[700]) for k,_ in selected) and len(selected)==6
