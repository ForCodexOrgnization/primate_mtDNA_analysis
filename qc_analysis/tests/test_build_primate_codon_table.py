import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(BIOPYTHON_AVAILABLE, 'Biopython is required for GenBank fixture tests')
class BuildPrimateCodonTableTests(unittest.TestCase):
    def make_record(self, path):
        # Plus CDS at 1..6 is ATGAAA. Minus CDS positions 7..12 encode ATGCCC.
        record = SeqRecord(Seq('ATGAAAGGGCAT'), id='TEST.1', name='TEST')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [
            SeqFeature(FeatureLocation(0, 6, strand=1), type='CDS', qualifiers={'gene':['ND1'], 'codon_start':['1'], 'transl_table':['2']}),
            SeqFeature(CompoundLocation([FeatureLocation(6, 9, strand=-1), FeatureLocation(9, 12, strand=-1)]), type='CDS', qualifiers={'gene':['COI'], 'codon_start':['1']}),
        ]
        SeqIO.write(record, path, 'genbank')

    def test_local_genbank_builds_coding_orientation_and_duplicates_samples(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); gbdir = d / 'gb'; gbdir.mkdir(); self.make_record(gbdir / 'TEST.1.gb')
            refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\taccession\tfamily\nS1\tSpecies one\tTEST.1\tFam\nS2\tSpecies one\tTEST.1\tFam\n')
            config = d / 'config.yaml'; output = d / 'table.tsv'; failures = d / 'failed.tsv'; summary = d / 'summary.tsv'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {gbdir}
    output_table: {output}
    failed_downloads_table: {failures}
    summary_table: {summary}
  settings:
    accession_columns: accession,reference_id
    skip_existing_genbank: true
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config)], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = list(csv.DictReader(output.open(), delimiter='\t'))
            self.assertEqual({r['sample'] for r in rows}, {'S1', 'S2'})
            s1 = [r for r in rows if r['sample'] == 'S1']
            self.assertEqual(len(s1), 12)
            plus = next(r for r in s1 if r['pos'] == '1')
            self.assertEqual((plus['codon_seq'], plus['codon_pos_in_triplet'], plus['codon_pos1_genomic']), ('ATG', '1', '1'))
            minus = next(r for r in s1 if r['pos'] == '12')
            self.assertEqual((minus['codon_seq'], minus['codon_pos_in_triplet']), ('ATG', '1'))
            self.assertEqual((minus['codon_pos1_genomic'], minus['codon_pos2_genomic'], minus['codon_pos3_genomic']), ('12', '11', '10'))
            self.assertEqual(minus['gene'], 'MT-CO1')

    def test_missing_accession_is_reported_without_stopping(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\nmissing\tSpecies\n')
            config = d / 'config.yaml'; failures = d / 'failed.tsv'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {d / 'gb'}
    output_table: {d / 'table.tsv'}
    failed_downloads_table: {failures}
    summary_table: {d / 'summary.tsv'}
  settings: {{accession_columns: accession, sample_column: sample, species_column: species}}
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config)], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            failure = next(csv.DictReader(failures.open(), delimiter='\t'))
            self.assertEqual(failure['sample'], 'missing')
            self.assertIn('No accession', failure['reason'])

    def test_mitos2_reference_fallback_selects_one_group_and_normalizes_rows(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            coordinate_fasta = d / 'Species_one.fa'; coordinate_fasta.write_text('>chrM\nATG\n')
            refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\taccession\nS1\tSpecies one\tACC.1\n')
            reference_rows = d / 'mitos_reference.tsv'
            fields = ['coordinate_reference_fasta', 'coordinate_reference_accession', 'gene', 'pos',
                      'codon_index', 'codon_pos_in_triplet', 'codon_pos1_genomic',
                      'codon_pos2_genomic', 'codon_pos3_genomic', 'codon_seq']
            first_group = [
                [coordinate_fasta, 'ACC.1', 'ND1', pos, '1', phase, '1', '2', '3', 'ATG']
                for pos, phase in [('1', '1'), ('2', '2'), ('3', '3')]
            ]
            # Duplicate the selected group and add a second accession-level group.
            second_group = [[d / 'other.fa', 'ACC.1', 'ND2', str(pos), '1', str(phase), '4', '5', '6', 'CCC']
                            for pos, phase in [(4, 1), (5, 2), (6, 3)]]
            with reference_rows.open('w', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t'); writer.writerow(fields)
                writer.writerows(first_group + first_group + second_group)
            output, summary, diagnostic = d / 'table.tsv', d / 'summary.tsv', d / 'fallback.tsv'
            config = d / 'config.yaml'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {d / 'gb'}
    species_fasta_dir: {d}
    output_table: {output}
    failed_downloads_table: {d / 'failed.tsv'}
    summary_table: {summary}
  settings:
    accession_columns: accession
    use_mitos2_if_genbank_fails: true
mitos2_annotation:
  paths:
    mitos2_reference_cds_table: {reference_rows}
    mitos2_fallback_selection_summary_table: {diagnostic}
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config)], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = list(csv.DictReader(output.open(), delimiter='\t'))
            self.assertEqual(len(rows), 3)
            self.assertEqual({row['gene'] for row in rows}, {'MT-ND1'})
            self.assertEqual({row['codon_pos_in_triplet'] for row in rows}, {'1', '2', '3'})
            summary_row = next(csv.DictReader(summary.open(), delimiter='\t'))
            self.assertEqual(summary_row['status'], 'completed_mitos2_fallback')
            self.assertIn('Built 3 coding-position rows for 1 samples (0 GenBank, 1 MITOS2 fallback; 0 failed, 0 other).', result.stdout)
            self.assertNotIn('Invalid codon_pos_in_triplet values detected.', summary_row['note'])
            self.assertIn('MITOS2 fallback duplicate rows collapsed: 3', summary_row['note'])
            diagnostic_row = next(csv.DictReader(diagnostic.open(), delimiter='\t'))
            self.assertEqual(diagnostic_row['fallback_match_mode'], 'coordinate_fasta')
            self.assertEqual(diagnostic_row['n_selected_rows_after_dedup'], '3')
            diagnostic_rows = list(csv.DictReader(diagnostic.open(), delimiter='\t'))
            self.assertEqual(len(diagnostic_rows), 2)
            rejected = next(row for row in diagnostic_rows if row['selection_status'] == 'rejected')
            self.assertEqual(rejected['coordinate_reference_fasta'], str(d / 'other.fa'))
            self.assertEqual(rejected['rejection_reason'], 'no_exact_coordinate_reference_fasta')
            self.assertEqual(list(csv.DictReader((d / 'failed.tsv').open(), delimiter='\t')), [])

    def test_mixed_genbank_mitos2_and_failure_is_deterministic(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); gbdir = d / 'gb'; gbdir.mkdir(); self.make_record(gbdir / 'TEST.1.gb')
            refs = d / 'refs.tsv'
            refs.write_text('sample\tspecies\taccession\nGB\tGenBank species\tTEST.1\nFB\tFallback species\nFAIL\tFailed species\n')
            mitos = d / 'mitos.tsv'
            mitos.write_text('sample\tspecies\tgene\tpos\tcodon_index\tcodon_pos_in_triplet\tcodon_pos1_genomic\tcodon_pos2_genomic\tcodon_pos3_genomic\tcodon_seq\n'
                             'FB\tFallback species\tND1\t1\t1\t1\t1\t2\t3\tATG\n')
            outputs = []
            for workers in (1, 4):
                output, summary, failures, diagnostic = (d / f'table-{workers}.tsv', d / f'summary-{workers}.tsv',
                                                          d / f'failed-{workers}.tsv', d / f'fallback-{workers}.tsv')
                config = d / f'config-{workers}.yaml'
                config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {gbdir}
    output_table: {output}
    failed_downloads_table: {failures}
    summary_table: {summary}
  settings:
    accession_columns: accession
    skip_existing_genbank: true
    use_mitos2_if_genbank_fails: true
mitos2_annotation:
  paths:
    mitos2_cds_table: {mitos}
    mitos2_fallback_selection_summary_table: {diagnostic}
''')
                result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config), '--workers', str(workers)], cwd=ROOT, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('Built 13 coding-position rows for 2 samples (1 GenBank, 1 MITOS2 fallback; 1 failed, 0 other).', result.stdout)
                outputs.append(tuple(path.read_text() for path in (output, summary, failures, diagnostic)))
                self.assertIn('  1 parse targets', result.stderr)
            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual({row['sample'] for row in csv.DictReader((d / 'table-1.tsv').open(), delimiter='\t')}, {'GB', 'FB'})
            self.assertEqual([row['sample'] for row in csv.DictReader((d / 'failed-1.tsv').open(), delimiter='\t')], ['FAIL'])

    def test_seven_mitos2_only_samples_need_no_genbank_parse_targets(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            samples = [f'S{i}' for i in range(7)]
            refs = d / 'refs.tsv'; refs.write_text('sample\tspecies\n' + ''.join(f'{sample}\tSpecies {sample}\n' for sample in samples))
            mitos = d / 'mitos.tsv'
            mitos.write_text('sample\tspecies\tgene\tpos\tcodon_index\tcodon_pos_in_triplet\tcodon_pos1_genomic\tcodon_pos2_genomic\tcodon_pos3_genomic\tcodon_seq\n' +
                             ''.join(f'{sample}\tSpecies {sample}\tND1\t1\t1\t1\t1\t2\t3\tATG\n' for sample in samples))
            config = d / 'config.yaml'; output, summary, failures = d / 'table.tsv', d / 'summary.tsv', d / 'failed.tsv'
            config.write_text(f'''build_primate_codon_table:
  paths:
    sample_ref_file: {refs}
    genbank_dir: {d / 'gb'}
    output_table: {output}
    failed_downloads_table: {failures}
    summary_table: {summary}
  settings:
    accession_columns: accession
    use_mitos2_if_genbank_fails: true
mitos2_annotation:
  paths:
    mitos2_cds_table: {mitos}
''')
            result = subprocess.run([sys.executable, str(ROOT / 'qc_analysis/scripts/build_primate_codon_table.py'), '--config', str(config), '--workers', '4'], cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('  0 parse targets', result.stderr)
            self.assertIn('Built 7 coding-position rows for 7 samples (0 GenBank, 7 MITOS2 fallback; 0 failed, 0 other).', result.stdout)
            self.assertEqual({row['sample'] for row in csv.DictReader(output.open(), delimiter='\t')}, set(samples))
            self.assertEqual({row['status'] for row in csv.DictReader(summary.open(), delimiter='\t')}, {'completed_mitos2_fallback'})
            self.assertEqual(list(csv.DictReader(failures.open(), delimiter='\t')), [])


class BuildPrimateCodonTableParallelHelperTests(unittest.TestCase):
    def setUp(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location('codon_builder', ROOT / 'qc_analysis/scripts/build_primate_codon_table.py')
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)

    def test_worker_resolution_uses_slurm_and_rejects_bad_values(self):
        from unittest.mock import patch
        with patch.dict('os.environ', {'SLURM_CPUS_PER_TASK': '7'}, clear=False):
            self.assertEqual(self.module.resolve_workers(None, {}, 'workers'), 7)
        with self.assertRaises(SystemExit) as error:
            self.module.resolve_workers('zero', {}, 'workers')
        self.assertIn('positive integer', str(error.exception))

    def test_fasta_index_matches_find_species_fasta(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); (d / 'Species_one.fa').write_text('>x\nA\n'); (d / 'Species_one.fasta.gz').write_text('not a real gzip')
            # The longer configured gzip extension takes precedence, as in the legacy lookup.
            index = self.module.build_fasta_index(d, '.fa,.fasta,.fasta.gz')
            self.assertEqual(self.module.find_species_fasta('Species one', d, '.fa,.fasta,.fasta.gz', index), self.module.find_species_fasta('Species one', d, '.fa,.fasta,.fasta.gz'))

    def test_atomic_failed_download_leaves_no_cache_file(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as td:
            destination = Path(td) / 'ACC.gb'
            class EmptyHandle:
                def __enter__(self): return self
                def __exit__(self, *args): pass
                def read(self): return ''
            import types
            self.module.Entrez = types.SimpleNamespace(efetch=None)
            with patch.object(self.module.Entrez, 'efetch', return_value=EmptyHandle()):
                with self.assertRaises(RuntimeError):
                    self.module.download('ACC', destination, {}, None)
            self.assertFalse(destination.exists())

    def test_summarize_build_status_counts_both_success_routes(self):
        counts = self.module.summarize_build_status([
            {'sample': 'genbank', 'status': 'completed'},
            {'sample': 'mitos2', 'status': 'completed_mitos2_fallback'},
            {'sample': 'failed', 'status': 'failed'},
        ])
        self.assertEqual(counts, {
            'total': 3, 'completed': 2, 'completed_genbank': 1,
            'completed_mitos2_fallback': 1, 'failed': 1, 'other': 0,
        })
        self.assertEqual(counts['completed'], counts['completed_genbank'] + counts['completed_mitos2_fallback'])

    def test_summarize_build_status_counts_seven_mitos2_successes(self):
        counts = self.module.summarize_build_status([
            {'sample': f'S{i}', 'status': 'completed_mitos2_fallback'} for i in range(7)
        ])
        self.assertEqual(counts['completed'], 7)
        self.assertEqual(counts['completed_genbank'], 0)
        self.assertEqual(counts['completed_mitos2_fallback'], 7)
        self.assertEqual(counts['failed'], 0)

    def test_fallback_group_selection_uses_all_priorities_and_marks_pre_lexical_ties(self):
        def rows(fasta, accession, genes, count):
            return [{'coordinate_reference_fasta': fasta, 'coordinate_reference_accession': accession,
                     'gene': gene, 'pos': str(index)}
                    for index, gene in enumerate((genes * ((count + len(genes) - 1) // len(genes)))[:count], 1)]

        genes = sorted(self.module.PROTEIN_CODING_GENES)
        candidates = (
            rows('/wrong.fa', 'ACC.1', genes, 11400) +
            rows('/sample.fa', 'WRONG', genes, 11400) +
            rows('/sample.fa', 'ACC.1', genes, 11400)
        )
        selected, mode, profiles, ambiguous = self.module.select_reference_fallback(
            candidates, '/sample.fa', 'ACC.1')

        self.assertEqual(mode, 'coordinate_fasta')
        self.assertEqual({row['coordinate_reference_accession'] for row in selected}, {'ACC.1'})
        self.assertFalse(ambiguous)
        self.assertEqual(len(profiles), 3)
        self.assertEqual(profiles[1]['rejection_reason'], 'no_exact_coordinate_reference_accession')

        tied = rows('/b.fa', 'B', genes, 11400) + rows('/a.fa', 'A', genes, 11400)
        selected, _, profiles, ambiguous = self.module.select_reference_fallback(tied, '', '')
        self.assertTrue(ambiguous)
        self.assertEqual(selected[0]['coordinate_reference_fasta'], '/a.fa')
        self.assertEqual(profiles[1]['rejection_reason'], 'deterministic_lexical_tiebreaker')

    def test_output_summary_consistency_warning(self):
        from unittest.mock import patch
        matching = [{'sample': 'S1', 'status': 'completed'}]
        self.assertTrue(self.module.warn_if_output_summary_disagree(matching, [{'sample': 'S1'}]))
        with patch('sys.stderr') as stderr:
            self.assertFalse(self.module.warn_if_output_summary_disagree(matching, [{'sample': 'S2'}]))
        message = ''.join(str(call.args[0]) for call in stderr.write.call_args_list)
        self.assertIn('WARNING: build summary and final codon table sample sets disagree.', message)
        self.assertIn('Successful-without-output: S1', message)
        self.assertIn('Output-without-success-status: S2', message)


if __name__ == '__main__':
    unittest.main()
