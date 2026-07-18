import csv, tempfile, unittest
from pathlib import Path
from qc_analysis.lib.mt_anchor_utils import rotate_sequence, sequence_sha256
from qc_analysis.scripts.discover_global_liftover_anchor import main as discover_main


def write_fa(p, name, seq): p.write_text(f'>{name}\n{seq}\n')
def rows(p):
    with p.open() as h: return list(csv.DictReader(h, delimiter='\t'))
def summary(p):
    return dict(line.split('\t', 1) for line in p.read_text().splitlines() if line)

class GlobalAnchorDiscoveryTests(unittest.TestCase):
    def run_discovery(self, manifest_rows, refs, min_len=20, header=True, species_dir=False, expect_error=None, human_seq=None):
        td=tempfile.TemporaryDirectory(); d=Path(td.name)
        human_seq=human_seq or 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'
        write_fa(d/'human.fa','chrM',human_seq)
        refdir=d/'refs'; refdir.mkdir()
        fasta_paths={}
        for name,seq in refs.items():
            fp=(refdir if species_dir else d)/f'{name}.fa'; write_fa(fp,'chrM',seq); fasta_paths[name]=fp
        man=d/'samples.tsv'
        with man.open('w') as h:
            if header:
                columns=manifest_rows[0] if manifest_rows and isinstance(manifest_rows[0], list) else ['sample','species','reference_id','species_fasta']
                if manifest_rows and isinstance(manifest_rows[0], list): manifest_rows=manifest_rows[1:]
                h.write('\t'.join(columns)+'\n')
                for row in manifest_rows:
                    vals=[]
                    for c in columns:
                        v=row.get(c,'') if isinstance(row, dict) else ''
                        if c == 'species_fasta' and v in fasta_paths: v=str(fasta_paths[v])
                        vals.append(str(v))
                    h.write('\t'.join(vals)+'\n')
            else:
                for sample,species in manifest_rows:
                    h.write(f'{sample}\t{species}\n')
        cfg=d/'cfg.yaml'; out=d/'ga'
        cfg.write_text(f'''global_anchor_discovery:\n  output_dir: {out}\n  sample_ref_file: {man}\n  human_fasta: {d/'human.fa'}\n  min_reference_length: {min_len}\n  max_reference_length: 100\n  max_n_fraction: 0.2\n  candidate_window_size: 5\n  min_anchor_column_occupancy: 0.66\n  min_window_mean_occupancy: 0.66\n  max_window_gap_fraction: 0.34\n  max_homopolymer_run: 20\n  aligner: mafft\n  aligner_options: "--auto --quiet"\n  allow_simple_alignment_fallback: true\ncoordinate_liftover:\n  paths:\n    species_fasta_dir: {refdir}\n    species_fasta_extensions: .fa,.fasta,.fna\n''')
        if expect_error:
            with self.assertRaisesRegex(RuntimeError, expect_error):
                discover_main(['--config', str(cfg)])
        else:
            discover_main(['--config', str(cfg)])
        return td,out,refdir

    def test_identical_circular_sequence_different_starts(self):
        base='AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'
        refs={'r1':rotate_sequence(base,1),'r2':rotate_sequence(base,9)}
        td,out,_=self.run_discovery([{'sample':'A','species':'sp','reference_id':'r1','species_fasta':'r1'},{'sample':'B','species':'sp','reference_id':'r2','species_fasta':'r2'}], refs)
        pos=rows(out/'reference_anchor_positions.tsv')
        self.assertEqual(len(pos),2)
        self.assertNotEqual(pos[0]['anchor_original_position'], pos[1]['anchor_original_position'])
        starts={rotate_sequence(refs[r['reference_id']], int(r['anchor_original_position']))[:8] for r in pos}
        self.assertEqual(len(starts),1)
        td.cleanup()

    def test_duplicate_reference_sequence_deduplicates(self):
        refs={'r1':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'}
        td,out,_=self.run_discovery([{'sample':'A','species':'sp','reference_id':'r1','species_fasta':'r1'},{'sample':'B','species':'sp','reference_id':'r1','species_fasta':'r1'}], refs)
        man=rows(out/'unique_reference_manifest.tsv')
        self.assertEqual(len(man),1); self.assertEqual(man[0]['sample_count'],'2')
        self.assertEqual(man[0]['sequence_sha256'], sequence_sha256(refs['r1']))
        td.cleanup()

    def test_low_quality_truncated_excluded_and_no_eligible_rejected(self):
        refs={'short':'ACGTACGT'}
        td,out,_=self.run_discovery([{'sample':'A','species':'short','reference_id':'short','species_fasta':'short'}], refs, min_len=20, expect_error='No eligible species references remained after QC')
        exc=rows(out/'excluded_references.tsv')
        self.assertIn('REFERENCE_TOO_SHORT', exc[0]['exclusion_reason'])
        td.cleanup()

    def test_deterministic_input_ordering(self):
        refs={'a':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT','b':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA'}
        td1,out1,_=self.run_discovery([{'sample':'A','species':'sp','reference_id':'a','species_fasta':'a'},{'sample':'B','species':'sp','reference_id':'b','species_fasta':'b'}], refs)
        td2,out2,_=self.run_discovery([{'sample':'B','species':'sp','reference_id':'b','species_fasta':'b'},{'sample':'A','species':'sp','reference_id':'a','species_fasta':'a'}], refs)
        self.assertEqual((out1/'global_anchor_selection.tsv').read_text(), (out2/'global_anchor_selection.tsv').read_text())
        td1.cleanup(); td2.cleanup()

    def test_headerless_two_column_manifest_resolves_species_dir(self):
        refs={'pan_troglodytes':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA'}
        td,out,_=self.run_discovery([('sample1','pan_troglodytes')], refs, header=False, species_dir=True)
        self.assertGreaterEqual(len(rows(out/'reference_anchor_positions.tsv')),1)
        s=summary(out/'global_anchor_summary.tsv')
        self.assertEqual(s['manifest_rows'],'1')
        self.assertEqual(s['resolved_sample_rows'],'1')
        td.cleanup()

    def test_explicit_species_fasta_still_supported(self):
        refs={'explicit':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA'}
        td,out,_=self.run_discovery([{'sample':'A','species':'does_not_match_dir','reference_id':'explicit','species_fasta':'explicit'}], refs)
        self.assertEqual(rows(out/'reference_anchor_positions.tsv')[0]['reference_id'],'explicit')
        td.cleanup()

    def test_missing_species_fasta_recorded_in_exclusions(self):
        td,out,_=self.run_discovery([('A','missing_species')], {}, header=False, species_dir=True, expect_error='No species references were resolved from sample_ref_file')
        exc=rows(out/'excluded_references.tsv')
        self.assertEqual(exc[0]['sample_names'],'A')
        self.assertIn('No species FASTA found', exc[0]['exclusion_reason'])
        td.cleanup()

    def test_human_only_msa_rejected(self):
        td,out,_=self.run_discovery([('A','missing_species')], {}, header=False, species_dir=True, expect_error='No species references were resolved from sample_ref_file')
        self.assertFalse((out/'reference_anchor_positions.tsv').exists())
        td.cleanup()

    def test_successful_output_contains_data_row_and_summary_counts(self):
        refs={'ok':'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTA'}
        td,out,_=self.run_discovery([{'sample':'A','species':'ok','reference_id':'ok','species_fasta':'ok'}], refs)
        self.assertGreaterEqual(len(rows(out/'reference_anchor_positions.tsv')),1)
        s=summary(out/'global_anchor_summary.tsv')
        self.assertEqual(s['unique_references'],'1')
        self.assertEqual(s['eligible_references'],'1')
        self.assertEqual(s['msa_sequence_count'],'2')
        self.assertEqual(s['reference_anchors_written'],'1')
        td.cleanup()

if __name__=='__main__': unittest.main()
