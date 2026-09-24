from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from qc_analysis.scripts.run_interspecies_contamination import (
    provenance_relationship,
    score_interpretation,
    specificity_weight,
)


def test_cross_species_specificity_weights():
    settings = {}
    assert specificity_weight(0.05, settings) == 1.0
    assert specificity_weight(0.08, settings) == 0.75
    assert specificity_weight(0.20, settings) == 0.50
    assert specificity_weight(0.40, settings) == 0.25
    assert specificity_weight(0.50, settings) == 0.0


def test_provenance_only_modifies_source_plausibility():
    settings = {
        "provenance_same_cohort_factor": 1.00,
        "provenance_same_project_factor": 0.90,
        "provenance_different_cohort_factor": 0.75,
        "provenance_different_project_factor": 0.60,
        "provenance_unknown_factor": 1.00,
    }
    provenance = {
        "A": {"project": "P1", "cohort": "C1"},
        "B": {"project": "P1", "cohort": "C1"},
        "C": {"project": "P1", "cohort": "C2"},
        "D": {"project": "P2", "cohort": "C3"},
    }
    assert provenance_relationship("A", "B", provenance, settings)["provenance_source_factor"] == 1.00
    assert provenance_relationship("A", "C", provenance, settings)["provenance_source_factor"] == 0.75
    assert provenance_relationship("A", "D", provenance, settings)["provenance_source_factor"] == 0.60
    assert provenance_relationship("A", "MISSING", provenance, settings)["provenance_source_factor"] == 1.00


def test_score_interpretation_bins():
    assert score_interpretation(None) == "not_scored_insufficient_evidence"
    assert score_interpretation(0.29) == "little_evidence"
    assert score_interpretation(0.30) == "weak_ambiguous_evidence"
    assert score_interpretation(0.50) == "candidate_evidence"
    assert score_interpretation(0.70) == "strong_evidence"


def test_production_wrapper_wires_full_native_artifact_workflow():
    wrapper = (ROOT / "qc_analysis/scripts/run_qc_preprocessing.sh").read_text()
    dedup = (ROOT / "qc_analysis/scripts/run_qc_preprocessing_deduplicated.sh").read_text()

    assert 'LOCAL_HET_QC_SCRIPT="qc_analysis/scripts/run_local_heteroplasmy_qc_with_expansion.sh"' in wrapper
    assert 'FINAL_FILTER_SCRIPT="qc_analysis/scripts/run_final_filter_with_heteroplasmy.py"' in wrapper
    assert "propagate_species_numt_positions.py" not in dedup
