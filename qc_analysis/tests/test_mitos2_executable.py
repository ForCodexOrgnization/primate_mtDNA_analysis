import importlib.util
from pathlib import Path
import pytest


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


def test_reference_tasks_are_one_based_per_target_species_and_report_completion(tmp_path):
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

    assert [row['task_id'] for row in rows] == [1, 2, 3]
    assert [row['reference_key'] for row in rows] == ['Species_one', 'Species_three', 'Species_two']
    assert [row['n_samples_using_reference'] for row in rows] == [1, 1, 1]
    assert [row['status'] for row in rows] == ['pending', 'pending', 'pending']


def test_references_prefers_target_variant_calling_fasta_for_cross_species_reference(tmp_path):
    module = load_module()
    fasta_dir = tmp_path / 'Ref_chrM'; fasta_dir.mkdir()
    target_fasta = fasta_dir / 'Target_species.fa'
    target_fasta.write_text('>chrM\nATG\n')
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(
        'target_species\tfinal_chrM_species\tfinal_chrM_accession\tfinal_chrM_refseq_accn\tchrM_expected_output_fasta\n'
        'Target_species\tOther species\tGB_1.1\tNC_123.4\t/references/chrM/independent/NC_123.4.fa\n'
    )
    samples = tmp_path / 'samples.tsv'; samples.write_text('sample\tspecies\nS1\tTarget_species\n')
    refs = module.references({'reference_manifest': str(manifest), 'sample_ref_file': str(samples),
                              'fasta_dir': str(fasta_dir), 'mitos2_raw_dir': str(tmp_path / 'raw')})

    ref, linked = refs[0]
    assert ref['mitos2_input_fasta'] == str(target_fasta)
    assert ref['coordinate_reference_fasta'] == str(target_fasta)
    assert ref['coordinate_reference_fasta_from_manifest'].endswith('NC_123.4.fa')
    assert ref['target_species'] == 'Target_species'
    assert ref['final_chrM_species'] == 'Other species'
    assert ref['final_chrM_accession'] == 'GB_1.1'
    assert ref['coordinate_reference_accession'] == 'GB_1.1'
    assert linked == [{'sample': 'S1', 'species': 'Target_species'}]


def test_references_sanitizes_manifest_fallback_and_skips_known_no_chrm(tmp_path):
    module = load_module()
    fallback = tmp_path / 'source.fa'; fallback.write_text('>NC_123.4 source\nATG\n')
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(
        'target_species\tfinal_reference_strategy\tchrM_expected_output_fasta\n'
        f'Fallback species\t\t{fallback}\n'
        'No chrM species\twg_only_no_chrM\t\n'
    )
    samples = tmp_path / 'samples.tsv'; samples.write_text('sample\tspecies\nS1\tFallback species\nS2\tNo chrM species\n')
    paths = {'reference_manifest': str(manifest), 'sample_ref_file': str(samples),
             'fasta_dir': str(tmp_path / 'Ref_chrM'), 'mitos2_raw_dir': str(tmp_path / 'raw')}
    refs = module.references(paths)
    fallback_ref = next(ref for ref, _ in refs if ref['target_species'] == 'Fallback species')
    skipped_ref = next(ref for ref, _ in refs if ref['target_species'] == 'No chrM species')

    assert Path(fallback_ref['mitos2_input_fasta']).read_text() == '>chrM\nATG\n'
    assert skipped_ref['initial_status'] == 'skipped_no_chrM_reference'
    assert next(row for row in module.task_rows(refs, paths) if row['target_species'] == 'No chrM species')['status'] == 'skipped_no_chrM_reference'


def test_missing_chrm_reference_is_skipped_not_fallback_task_and_is_summarized(tmp_path):
    module = load_module()
    manifest = tmp_path / 'manifest.tsv'
    samples = tmp_path / 'samples.tsv'
    manifest.write_text(
        'target_species\tfinal_chrM_species\tfinal_chrM_accession\tfinal_chrM_refseq_accn\tfinal_chrM_genbank_accn\tchrM_expected_output_fasta\tchrM_selection_status\tfinal_reference_strategy\treference_pairing_status\n'
        'WG only species\t\t\t\t\t\tmissing_chrM_ref\twg_only_no_chrM\tno_chrM_pair\n'
    )
    samples.write_text('sample\tspecies\nS1\tWG only species\n')
    paths = {
        'reference_manifest': str(manifest),
        'sample_ref_file': str(samples),
        'fasta_dir': str(tmp_path / 'fallback'),
        'mitos2_raw_dir': str(tmp_path / 'raw'),
    }

    refs = module.references(paths)

    assert len(refs) == 1
    ref, linked = refs[0]
    assert ref['status'] == 'skipped_no_chrM_reference'
    assert ref['coordinate_reference_fasta'] == ''
    assert module.task_rows(refs, paths) == []
    summary = module.collect_reference(ref, linked, paths, {})['summary_row']
    assert summary['status'] == 'skipped_no_chrM_reference'
    assert summary['chrM_selection_status'] == 'missing_chrM_ref'
    assert summary['final_reference_strategy'] == 'wg_only_no_chrM'
    assert summary['reference_pairing_status'] == 'no_chrM_pair'


def test_collect_reference_separates_reference_and_sample_codon_counts(tmp_path, monkeypatch):
    module = load_module()
    fasta = tmp_path / 'NC_002764.1.fa'
    fasta.write_text('>NC_002764.1\n' + 'ATG' * 10 + '\n')
    raw = tmp_path / 'raw' / 'NC_002764.1'
    raw.mkdir(parents=True)
    raw.joinpath('result.gff').write_text(
        'chrM\tmitos\tgene\t1\t30\t.\t+\t.\tID=gene_nad1;Name=nad1;gene_id=nad1\n'
    )
    ref = {
        'reference_key': 'NC_002764.1',
        'reference_species': 'Macaca_sylvanus',
        'coordinate_reference_accession': 'NC_002764.1',
        'coordinate_reference_fasta': str(fasta),
    }
    paths = {'mitos2_raw_dir': str(tmp_path / 'raw')}
    monkeypatch.setattr(module, 'build_reference_codon_rows', lambda *args: [{'pos': i} for i in range(30)])

    result = module.collect_reference(ref, [{'sample': 'sample-1', 'species': 'Macaca sylvanus'}], paths, {'genetic_code': 2})

    assert result['status'] == 'completed'
    assert len(result['features']) == 1
    assert len(result['reference_codon_rows']) == 30
    assert len(result['sample_codon_rows']) == 30
    summary = result['summary_row']
    assert summary['n_cds_features'] == 1
    assert summary['n_linked_samples'] == 1
    assert summary['n_reference_coding_position_rows'] == 30
    assert summary['n_sample_level_coding_position_rows'] == 30
    assert summary['result_gff_exists'] is True


def test_collect_reference_reports_parser_cds_detection_failure(tmp_path):
    module = load_module()
    raw = tmp_path / 'raw' / 'NC_002764.1'
    raw.mkdir(parents=True)
    # The diagnostics recognizes nad1 as coding even when parsing produces no features.
    raw.joinpath('result.gff').write_text('chrM\tmitos\tgene\t1\t30\t.\t+\t.\tName=nad1\n')
    ref = {'reference_key': 'NC_002764.1', 'reference_species': 'Macaca_sylvanus',
           'coordinate_reference_accession': 'NC_002764.1', 'coordinate_reference_fasta': str(tmp_path / 'missing.fa')}

    # Make the parser intentionally return no rows to exercise the diagnostic branch.
    module.parse_outputs = lambda raw_dir, reference: ([], [])
    result = module.collect_reference(ref, [], {'mitos2_raw_dir': str(tmp_path / 'raw')}, {})

    assert result['status'] == 'failed_parser_cds_gene_detection'
    assert result['summary_row']['n_cds_features'] == 0
    assert result['summary_row']['n_gff_cds_like_gene_rows'] == 1


def test_reference_codon_rows_select_gff_seqid_and_isolate_bad_gene(tmp_path):
    module = load_module()
    if module.SeqIO is None:
        pytest.skip('Biopython is required for FASTA codon-row coverage')
    fasta = tmp_path / 'multi.fa'
    fasta.write_text('>wrong\n' + 'A' * 30 + '\n>chrM\n' + 'ATG' * 10 + '\n')
    raw = tmp_path / 'raw'; raw.mkdir()
    ref = {'reference_key': 'NC_002764.1', 'coordinate_reference_accession': 'NC_002764.1', 'raw_dir': str(raw)}
    features = [
        {'feature_type': 'CDS', 'gff_seqid': 'chrM', 'start': '1', 'end': '30', 'strand': '+', 'gene': 'MT-ND1', 'gene_raw': 'nad1'},
        {'feature_type': 'CDS', 'gff_seqid': 'chrM', 'start': '31', 'end': '40', 'strand': '+', 'gene': 'MT-ND2', 'gene_raw': 'nad2'},
    ]
    rows = module.build_reference_codon_rows(features, fasta, ref, '2')
    debug = list(__import__('csv').DictReader((raw / 'mitos2_reference_codon_debug.tsv').open(), delimiter='\t'))
    assert len(rows) == 30
    assert rows[0]['seq_name'] == 'chrM'
    assert [row['status'] for row in debug] == ['completed', 'failed']
    assert debug[1]['error']
