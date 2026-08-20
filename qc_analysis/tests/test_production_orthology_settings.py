from pathlib import Path

from qc_analysis.lib.simple_yaml import read_simple_yaml


ROOT=Path(__file__).resolve().parents[2]


def test_production_trna_and_rrna_orthology_settings():
    config=read_simple_yaml(ROOT/'config/qc_preprocessing.yaml')
    trna=config['trna_match']['settings']
    rrna=config['rrna_match']['settings']
    assert trna['require_compensated_for_strict_stem'] is False
    assert trna['strict_stem_require_reference_pair_type_match'] is False
    assert rrna['require_same_rrna_element'] is False
    assert rrna['high_conf_loop_max_frac_delta'] == 0.002
