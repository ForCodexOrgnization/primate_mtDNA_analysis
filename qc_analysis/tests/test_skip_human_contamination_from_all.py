from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "qc_analysis" / "scripts" / "run_qc_preprocessing.sh"


def wrapper_text() -> str:
    return WRAPPER.read_text(encoding="utf-8")


def section(text: str, start: str, end: str) -> str:
    left = text.index(start)
    right = text.index(end, left)
    return text[left:right]


def test_submit_all_excludes_human_contamination_but_keeps_standalone_step():
    text = wrapper_text()
    submit = section(text, "submit_workflow() {", "if [[ \"$ARRAY_TASK_MODE\" == \"1\" ]]")

    assert "human_contamination" not in submit
    assert "human_contamination)" in text


def test_direct_all_excludes_human_contamination():
    text = wrapper_text()
    direct_all = section(text, "  all)\n", "    ;;\nesac")

    assert '"$HUMAN_CONTAMINATION_SCRIPT"' not in direct_all
    assert '"$PRIMATE_BACKGROUND_SCRIPT"' in direct_all
    assert '"$FINAL_FILTER_SCRIPT"' in direct_all


def test_help_documents_human_contamination_as_standalone_only():
    text = wrapper_text()

    assert "human_contamination              Screen annotated human-coordinate alleles with primate background correction (standalone; temporarily excluded from all)." in text
    assert "all                              Run production preprocessing/annotation steps, currently excluding human_contamination." in text
