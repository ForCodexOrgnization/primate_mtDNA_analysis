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
            self.assertNotIn('Invalid codon_pos_in_triplet values detected.', summary_row['note'])
            self.assertIn('MITOS2 fallback duplicate rows collapsed: 3', summary_row['note'])
            diagnostic_row = next(csv.DictReader(diagnostic.open(), delimiter='\t'))
            self.assertEqual(diagnostic_row['fallback_match_mode'], 'coordinate_fasta')
            self.assertEqual(diagnostic_row['n_selected_rows_after_dedup'], '3')


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


if __name__ == '__main__':
    unittest.main()
