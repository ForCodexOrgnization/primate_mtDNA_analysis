from pathlib import Path
import pytest
from qc_analysis.lib.simple_yaml import read_simple_yaml

def test_nested_scalars_quotes_comments_and_nulls(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text('''intraspecies_contamination:\n  enabled: true\n  vcf_dir: null\n  outdir: "results/#keep" # comment\n  thresholds: {}\n  dp_min: 100\n  low: 0.01\n''')
    value = read_simple_yaml(path)["intraspecies_contamination"]
    assert value == {"enabled": True, "vcf_dir": None, "outdir": "results/#keep", "thresholds": {}, "dp_min": 100, "low": .01}

def test_malformed_indentation_fails(tmp_path):
    path = tmp_path / "bad.yaml"; path.write_text("value: true\n  child: false\n")
    with pytest.raises(ValueError, match="indentation"):
        read_simple_yaml(path)
