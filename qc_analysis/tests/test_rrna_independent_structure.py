import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from qc_analysis.lib.match_utils import info_parse, yaml
from qc_analysis.scripts.run_rrna_match import element_match, match_tier


ROOT = Path(__file__).resolve().parents[2]
PRODUCTION_LOOP_THRESHOLD = float(
    yaml(ROOT / "config/qc_preprocessing.yaml")["rrna_match"]["settings"][
        "high_conf_loop_max_frac_delta"
    ]
)


def record_info(path):
    record = next(line for line in path.read_text().splitlines() if not line.startswith("#"))
    return info_parse(record.split("\t")[7])


def write_common_inputs(d, human_rows, species_rows, partner_human_pos="10"):
    human_regions = d / "human_regions.tsv"
    species_regions = d / "species_regions.tsv"
    human_structure = d / "human_structure.tsv"
    species_structure = d / "species_structure.tsv"
    sample_map = d / "sample_reference_map.tsv"
    maps = d / "maps"
    maps.mkdir()

    human_regions.write_text("chrom\tstart\tend\trrna_gene\tstrand\nchrM\t1\t20\tMT-RNR1\t+\n")
    species_regions.write_text(
        "reference_key\treference_species\tcoordinate_reference_sequence_sha256\tstart\tend\trrna_gene\tstrand\n"
        "wrong_ref\tSpecies one\twrong_sha\t1\t20\tMT-RNR2\t+\n"
        "ref1\tSpecies one\tsha1\t1\t20\tMT-RNR1\t+\n"
    )
    human_structure.write_text(
        "rrna_gene\thuman_pos\tlocal_pos\tbase\tstruct_class\tstruct_element\t"
        "paired_human_pos\tpaired_local_pos\tpaired_base\tpair_type\tpair_state\n"
        + "\n".join(human_rows) + "\n"
    )
    species_structure.write_text(
        "reference_key\treference_species\tcoordinate_reference_accession\t"
        "coordinate_reference_fasta\tcoordinate_reference_sequence_sha256\t"
        "rrna_gene\tgenomic_pos\tlocal_pos\tbase\tstruct_class\tstruct_element\tpaired_genomic_pos\t"
        "paired_local_pos\tpaired_base\tpair_type\tpair_state\tannotation_source\tstructure_source\n"
        + "\n".join(species_rows) + "\n"
    )
    sample_map.write_text(
        "sample\tspecies\tspecies_key\treference_key\tcoordinate_reference_fasta\t"
        "coordinate_reference_accession\tcoordinate_reference_sequence_sha256\n"
        "S1\tSpecies one\tspecies_one\tref1\t/ref.fa\tACC.1\tsha1\n"
    )
    (maps / "S1.coordinate_map.tsv").write_text(
        "species_pos_original\thuman_pos_canonical\n"
        f"1\t1\n11\t{partner_human_pos}\n12\t11\n"
    )
    input_vcf = d / "input.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chrM\t1\t.\tG\tA\t.\tPASS\tSRC_CHROM=species;SRC_POS=1;SRC_REF=A;SRC_ALT=A;MTLIFT_HUMAN_POS=1\n"
    )
    config = d / "config.yaml"
    config.write_text(f"""rrna_match:
  paths:
    input_vcf_dir: {d}
    fallback_codon_vcf_dir: {d}
    fallback_raw_vcf_dir: {d}
    output_dir: {d}
    reports_dir: {d / 'reports'}
    coordinate_map_dir: {maps}
    human_rrna_table: {human_regions}
    species_rrna_table: {species_regions}
    human_rrna_structure_table: {human_structure}
    species_rrna_structure_table: {species_structure}
    sample_reference_map: {sample_map}
  settings:
    input_vcf_pattern: "{{sample}}.vcf"
    fallback_codon_vcf_pattern: "{{sample}}.vcf"
    fallback_raw_vcf_pattern: "{{sample}}.vcf"
    output_suffix: ".out.vcf"
    use_rrna_structure_table: true
""")
    return config, input_vcf


def run_rrna_case(human_rows, species_rows, partner_human_pos="10"):
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        config, input_vcf = write_common_inputs(d, human_rows, species_rows, partner_human_pos)
        output = d / "out.vcf"
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "qc_analysis/scripts/run_rrna_match.py"),
                "--config", str(config),
                "--sample", "S1",
                "--input", str(input_vcf),
                "--output", str(output),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        assert completed.returncode == 0, completed.stderr
        return record_info(output)


def test_stem_stem_compensatory_pair_type_is_structurally_conserved():
    info = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tstem\tH1\t10\t10\tC\tG-C\tcanonical"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tstem\tH1\t11\t11\tT\tA-U\tcanonical\tMITOS2\ttoy.sto"],
    )

    assert info["MTRRNA_STRUCTURE_MATCH"] == "STEM_STEM"
    assert info["MTRRNA_PAIR_RELATION_MATCH"] == "yes"
    assert info["MTRRNA_MATCH_TIER"] == "HIGH_CONF_STEM"
    assert info["MTRRNA_H_PAIR_TYPE"] == "G-C"
    assert info["MTRRNA_S_PAIR_TYPE"] == "A-U"
    assert info["MTRRNA_H_ALT_EFFECT"] == "canonical_to_noncanonical"


def test_stem_to_loop_is_discordant_without_projecting_human_structure():
    info = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tstem\tH1\t10\t10\tC\tG-C\tcanonical"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tloop\tH1\t.\t.\t.\t.\tunpaired\tMITOS2\ttoy.sto"],
    )

    assert info["MTRRNA_STRUCTURE_MATCH"] == "STEM_LOOP"
    assert info["MTRRNA_PAIR_RELATION_MATCH"] == "NA"
    assert info["MTRRNA_MATCH_TIER"] == "STRUCTURE_DISCORDANT"


def test_loop_loop_is_high_conf_loop():
    info = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tloop\tH1\t.\t.\t.\t.\tunpaired"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tloop\tH1\t.\t.\t.\t.\tunpaired\tMITOS2\ttoy.sto"],
    )

    assert info["MTRRNA_STRUCTURE_MATCH"] == "LOOP_LOOP"
    assert info["MTRRNA_ELEMENT_MATCH"] == "yes"
    assert info["MTRRNA_MATCH_TIER"] == "HIGH_CONF_LOOP"


def test_missing_species_structure_is_unknown_not_human_projected():
    info = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tstem\tH1\t10\t10\tC\tG-C\tcanonical"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t2\t2\tA\tloop\tH1\t.\t.\t.\t.\tunpaired\tMITOS2\ttoy.sto"],
    )

    assert info["MTRRNA_H_CLASS"] == "stem"
    assert info["MTRRNA_S_CLASS"] == "unknown"
    assert info["MTRRNA_STRUCTURE_MATCH"] == "UNKNOWN"
    assert info["MTRRNA_MATCH_TIER"] == "STRUCTURE_UNKNOWN"


def test_wrong_species_partner_does_not_receive_high_conf_stem():
    info = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tstem\tH1\t10\t10\tC\tG-C\tcanonical"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tstem\tH1\t12\t12\tT\tA-U\tcanonical\tMITOS2\ttoy.sto"],
    )

    assert info["MTRRNA_STRUCTURE_MATCH"] == "STEM_STEM"
    assert info["MTRRNA_PAIR_RELATION_MATCH"] == "no"
    assert info["MTRRNA_MATCH_TIER"] != "HIGH_CONF_STEM"


@pytest.mark.parametrize(
    "structure,pair_relation,elements,delta,require_element,expected",
    [
        ("LOOP_LOOP", "NA", "yes", 0.002, True, "HIGH_CONF_LOOP"),
        ("LOOP_LOOP", "NA", "yes", 0.003, True, "LOW_CONF"),
        ("LOOP_LOOP", "NA", "no", 0.002, True, "LOW_CONF"),
        ("LOOP_LOOP", "NA", "no", 0.002, False, "HIGH_CONF_LOOP"),
        ("LOOP_LOOP", "NA", "yes", None, True, "LOW_CONF"),
        ("STEM_STEM", "yes", "yes", 0.99, True, "HIGH_CONF_STEM"),
        ("STEM_STEM", "yes", "no", 0.01, True, "LOW_CONF"),
        ("STEM_STEM", "yes", "no", 0.01, False, "HIGH_CONF_STEM"),
        ("STEM_STEM", "no", "yes", 0.01, True, "LOW_CONF"),
        ("STEM_LOOP", "NA", "yes", 0.01, True, "STRUCTURE_DISCORDANT"),
        ("LOOP_STEM", "NA", "yes", 0.01, True, "STRUCTURE_DISCORDANT"),
        ("UNKNOWN", "NA", "yes", 0.01, True, "STRUCTURE_UNKNOWN"),
    ],
)
def test_match_tier_uses_element_and_loop_fraction_requirements(
    structure, pair_relation, elements, delta, require_element, expected
):
    assert match_tier(
        "OK", True, structure, pair_relation, elements, delta, require_element,
        PRODUCTION_LOOP_THRESHOLD,
    ) == expected


@pytest.mark.parametrize(
    "structure,pair_relation,delta,expected",
    [
        ("STEM_STEM", "yes", 0.99, "HIGH_CONF_STEM"),
        ("STEM_STEM", "no", 0.01, "LOW_CONF"),
        ("LOOP_LOOP", "NA", 0.002, "HIGH_CONF_LOOP"),
        ("LOOP_LOOP", "NA", 0.003, "LOW_CONF"),
        ("LOOP_LOOP", "NA", None, "LOW_CONF"),
    ],
)
def test_default_tiers_do_not_require_unknown_elements(
    structure, pair_relation, delta, expected
):
    assert match_tier(
        "OK", True, structure, pair_relation, ".", delta, False,
        PRODUCTION_LOOP_THRESHOLD,
    ) == expected


@pytest.mark.parametrize(
    "delta,expected",
    [
        (0.002, "HIGH_CONF_LOOP"),
        (0.0020001, "LOW_CONF"),
        (0.0025, "LOW_CONF"),
        (0.003, "LOW_CONF"),
        (None, "LOW_CONF"),
    ],
)
def test_loop_loop_uses_production_fraction_threshold(delta, expected):
    assert PRODUCTION_LOOP_THRESHOLD == 0.002
    assert match_tier(
        "OK", True, "LOOP_LOOP", "NA", "yes", delta, False,
        PRODUCTION_LOOP_THRESHOLD,
    ) == expected


def test_loop_fraction_threshold_does_not_apply_to_high_conf_stem():
    assert match_tier(
        'OK',True,'STEM_STEM','yes','.',0.99,False,0.002
    ) == 'HIGH_CONF_STEM'


def test_unknown_elements_produce_high_conf_tiers_with_production_default():
    stem = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tstem\t.\t10\t10\tC\tG-C\tcanonical"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tstem\t.\t11\t11\tT\tA-U\tcanonical\tMITOS2\ttoy.sto"],
    )
    loop = run_rrna_case(
        ["MT-RNR1\t1\t1\tG\tloop\t.\t.\t.\t.\t.\tunpaired"],
        ["ref1\tSpecies\tACC\t/ref.fa\tsha1\tMT-RNR1\t1\t1\tA\tloop\t.\t.\t.\t.\t.\tunpaired\tMITOS2\ttoy.sto"],
    )
    assert stem["MTRRNA_ELEMENT_MATCH"] == loop["MTRRNA_ELEMENT_MATCH"] == "."
    assert stem["MTRRNA_MATCH_TIER"] == "HIGH_CONF_STEM"
    assert loop["MTRRNA_MATCH_TIER"] == "HIGH_CONF_LOOP"


def test_element_match_trims_only_and_preserves_missing_values():
    assert element_match(" H1 ", "H1") == "yes"
    assert element_match("H1", "H2") == "no"
    assert element_match("h1", "H1") == "no"
    assert element_match(".", "H1") == "."
