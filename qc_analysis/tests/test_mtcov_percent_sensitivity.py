from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "qc_analysis/scripts/run_mtcov_percent_sensitivity.py"


def load_module():
    spec = spec_from_file_location("mtcov_percent_sensitivity", SCRIPT)
    mod = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def test_six_matched_mtcov_percent_scenarios():
    mod = load_module()
    rows = mod.scenarios()
    assert len(rows) == 6
    observed = {
        (r["mtcov"], r["percent_depth"], r["fraction"], r["nuclear"], r["mtcn"])
        for r in rows
    }
    expected = {
        (60, 60, 85, 5, 40),
        (60, 60, 90, 5, 40),
        (80, 80, 85, 5, 40),
        (80, 80, 90, 5, 40),
        (100, 100, 85, 5, 40),
        (100, 100, 90, 5, 40),
    }
    assert observed == expected


def test_percent_at_depth_uses_greater_than_or_equal(tmp_path):
    mod = load_module()
    cov = tmp_path / "coverage.tsv"
    cov.write_text(
        "chrom\tpos\ttarget\tcoverage\n"
        "chrM\t1\tmt\t59\n"
        "chrM\t2\tmt\t60\n"
        "chrM\t3\tmt\t80\n"
        "chrM\t4\tmt\t100\n"
    )
    assert mod.percent_at_depth(cov, 60) == 75.0
    assert mod.percent_at_depth(cov, 80) == 50.0
    assert mod.percent_at_depth(cov, 100) == 25.0
