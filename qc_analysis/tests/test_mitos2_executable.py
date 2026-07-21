import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'qc_analysis/scripts/run_mitos2_annotation.py'


def load_module():
    spec = importlib.util.spec_from_file_location('run_mitos2_annotation', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_mitos2_uses_runmitos_after_conda_activation(monkeypatch):
    module = load_module()
    captured = {}

    class Result:
        returncode = 0
        stdout = 'CONDA_PREFIX=/tmp/mitos2\nMITOS executable=/tmp/mitos2/bin/runmitos\nUsing MITOS2 executable: /tmp/mitos2/bin/runmitos\n'
        stderr = ''

    def fake_run(args, **kwargs):
        captured['command'] = args[2]
        return Result()

    monkeypatch.setattr(module.subprocess, 'run', fake_run)
    executable, output = module.command({'conda_module': 'miniconda', 'conda_env': 'mitos2'})

    assert executable == 'runmitos'
    assert output == Result.stdout
    assert 'module load miniconda' in captured['command']
    assert 'source "$(conda info --base)/etc/profile.d/conda.sh"' in captured['command']
    assert 'conda activate mitos2' in captured['command']
    assert 'command -v runmitos' in captured['command']
    assert 'command -v ' + 'MITOS' + '2' not in captured['command']
    assert 'command -v ' + 'runmitos' + '.py' not in captured['command']


def test_mitos2_command_templates_invoke_runmitos():
    module = load_module()
    commands = module.templates('runmitos', 'input.fa', 'output', {'genetic_code': 2, 'threads': 4})

    assert all(command.startswith('runmitos ') for command in commands)
    assert commands == [
        'runmitos -i input.fa -c 2 -o output -r refseq81m --best --noplots',
        'runmitos --input input.fa -c 2 -o output -r refseq81m --best --noplots',
    ]
    assert '-I ' not in ' '.join(commands)
    assert '--fasta' not in ' '.join(commands)
    assert '--threads' not in ' '.join(commands)


def test_mitos2_command_templates_quote_paths_and_include_optional_refdir():
    module = load_module()
    commands = module.templates(
        Path('/opt/mitos/bin/runmitos'),
        Path('/tmp/input fasta.fa'),
        Path('/tmp/output directory'),
        {'genetic_code': 2, 'refseqver': 'refseq99m', 'refdir': '/tmp/reference data'},
    )

    assert commands[0] == (
        "/opt/mitos/bin/runmitos -i '/tmp/input fasta.fa' -c 2 "
        "-o '/tmp/output directory' -r refseq99m -R '/tmp/reference data' "
        "--best --noplots"
    )
    assert commands[1].startswith("/opt/mitos/bin/runmitos --input '/tmp/input fasta.fa'")


def test_parse_outputs_prefers_mitos_gff_and_normalizes_mitos_names(tmp_path):
    module = load_module()
    (tmp_path / 'result.gff').write_text(
        'chrM\tmitos\tregion\t1\t16965\t.\t+\t.\tID=chrM:1..16965\n'
        'chrM\tmitfi\tncRNA_gene\t1\t66\t.\t+\t.\tID=gene_trnF;Name=trnF;gene_id=trnF\n'
        'chrM\tmitfi\ttRNA\t1\t66\t.\t+\t.\tID=transcript_trnF(gaa);Name=trnF(gaa)\n'
        'chrM\tmitfi\trRNA\t67\t1022\t.\t+\t.\tID=transcript_rrnS;Name=rrnS\n'
        'chrM\tmitos\tgene\t2738\t3694\t.\t+\t.\tID=gene_nad1;Name=nad1;gene_id=nad1\n'
        'chrM\tmitfi\texon\t2738\t3694\t.\t+\t.\tParent=transcript_nad1;Name=nad1\n'
    )
    # This fallback must not be read while the GFF produces valid features.
    (tmp_path / 'result.bed').write_text('chrM\t1\t10\tnad2\t0\t+\n')
    features, diagnostics = module.parse_outputs(tmp_path, {'reference_key': 'test'})

    assert [feature['feature_type'] for feature in features] == ['tRNA', 'rRNA', 'CDS']
    assert [feature['gene'] for feature in features] == ['trnF', 'MT-RNR1', 'MT-ND1']
    assert features[-1]['gene_raw'] == 'nad1'
    assert len(diagnostics) == 1
    assert diagnostics[0]['file'].endswith('result.gff')
    assert module.gff_diagnostics(tmp_path) == {
        'result_gff_exists': True,
        'n_gff_gene_rows': 1,
        'n_gff_cds_like_gene_rows': 1,
        'n_gff_trna_rows': 1,
        'n_gff_rrna_rows': 1,
    }


def test_parser_failure_status_flags_unrecognized_cds_like_gff_genes():
    module = load_module()
    diagnostics = {
        'result_gff_exists': True,
        'n_gff_gene_rows': 1,
        'n_gff_cds_like_gene_rows': 1,
        'n_gff_trna_rows': 0,
        'n_gff_rrna_rows': 0,
    }

    assert module.parser_failure_status([], diagnostics) == 'failed_parser_cds_gene_detection'


def test_reference_tasks_are_one_based_deduplicated_and_report_completion(tmp_path):
    module = load_module()
    manifest = tmp_path / 'manifest.tsv'
    samples = tmp_path / 'samples.tsv'
    raw = tmp_path / 'raw'
    manifest.write_text(
        'target_species\tfinal_chrM_species\tfinal_chrM_accession\tchrM_expected_output_fasta\n'
        'Species one\tReference one\tACC.1\t/tmp/ref-one.fa\n'
        'Species two\tReference one\tACC.1\t/tmp/ref-one.fa\n'
        'Species three\tReference three\tACC.3\t/tmp/ref-three.fa\n'
    )
    samples.write_text('sample\tspecies\nS1\tSpecies one\nS2\tSpecies two\nS3\tSpecies three\n')
    paths = {'reference_manifest': str(manifest), 'sample_ref_file': str(samples), 'mitos2_raw_dir': str(raw)}
    refs = module.references(paths)
    (raw / 'ACC.1' / 'mitos2.completed.ok').parent.mkdir(parents=True)
    (raw / 'ACC.1' / 'mitos2.completed.ok').write_text('completed\n')

    rows = module.task_rows(refs, paths)

    assert [row['task_id'] for row in rows] == [1, 2]
    assert [row['reference_key'] for row in rows] == ['ACC.1', 'ACC.3']
    assert [row['n_samples_using_reference'] for row in rows] == [2, 1]
    assert [row['status'] for row in rows] == ['completed', 'pending']
