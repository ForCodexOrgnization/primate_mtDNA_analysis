import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from qc_analysis.lib.match_utils import info_parse
from qc_analysis.scripts.run_codon_match import mutate_codon


ROOT = Path(__file__).resolve().parents[2]


def record_info(path):
    record = next(line for line in path.read_text().splitlines() if not line.startswith('#'))
    return info_parse(record.split('\t')[7])


def run_script(script, config, sample, input_vcf, output):
    return subprocess.run(
        [sys.executable, str(ROOT / 'qc_analysis/scripts' / script), '--config', str(config),
         '--sample', sample, '--input', str(input_vcf), '--output', str(output)],
        cwd=ROOT, text=True, capture_output=True,
    )


class MatchScriptSmokeTests(unittest.TestCase):
    def write_vcf(self, directory, info, alt='C'):
        path = directory / 'input.vcf'
        path.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                        f'chrM\t10\t.\tT\t{alt}\t.\tPASS\t{info}\n')
        return path

    def test_codon_uses_source_alt_to_construct_alt_codon(self):
        self.assertEqual(mutate_codon('acg', 2, 'T'), 'ATG')
        self.assertEqual(mutate_codon('ACG', 4, 'T'), '.')
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); primate = d / 'primate.tsv'; human = d / 'human.tsv'
            primate.write_text('sample\tpos\tgene\tcodon_seq\tcodon_pos_in_triplet\nS1\t10\tMT-ND1\tACG\t2\n')
            human.write_text('pos\tgene\tcodon_seq\tcodon_pos_in_triplet\n10\tMT-ND1\tATG\t2\n')
            config = d / 'config.yaml'
            config.write_text(f'''codon_match:
  paths:
    input_vcf_dir: {d}
    output_dir: {d}
    reports_dir: {d / 'reports'}
    all_primate_position_codon_table: {primate}
    human_codon_table: {human}
  settings:
    strict_phase_match: true
    input_vcf_pattern: "{{sample}}.vcf"
    output_suffix: ".out.vcf"
''')
            output = d / 'codon.vcf'
            result = run_script('run_codon_match.py', config, 'S1', self.write_vcf(d, 'SRC_CHROM=species;SRC_POS=10;SRC_REF=C;SRC_ALT=T'), output)
            self.assertEqual(result.returncode, 0, result.stderr)
            info = record_info(output)
            self.assertEqual(info['MTCODON_PRIMATE_ALT_CODON'], 'ATG')
            self.assertEqual(info['MTCODON_STATUS'], 'PASS')

    def test_trna_region_match_compares_structural_classes_not_ids(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); human = d / 'human.tsv'; species = d / 'S1.tsv'
            header = 'chrom\tpos\ttrna_id\tlocal_pos\tstruct_class\tstruct_element\tpair_type\tpair_state\tpaired_local_pos\tpaired_genomic_pos\tpaired_base\n'
            human.write_text(header + 'chrM\t10\tTRNA-H\t1\tstem\tacceptor\tWC\tpaired\t2\t20\tG\n')
            species.write_text(header + 'species\t10\tTRNA-S\t1\tstem\tacceptor\tWC\tpaired\t2\t20\tG\n')
            config = d / 'config.yaml'
            config.write_text(f'''trna_match:
  paths:
    input_vcf_dir: {d}
    fallback_input_vcf_dir: {d}
    output_dir: {d}
    reports_dir: {d / 'reports'}
    coordinate_map_dir: {d / 'maps'}
    human_trna_index: {human}
    species_trna_index_dir: {d}
    species_trna_index_template: "{{species_trna_index_dir}}/{{sample}}.tsv"
  settings:
    input_vcf_pattern: "{{sample}}.vcf"
    fallback_input_vcf_pattern: "{{sample}}.vcf"
    output_suffix: ".out.vcf"
    species_trna_coord_space: original
    require_compensated_for_strict_stem: false
''')
            output = d / 'trna.vcf'
            result = run_script('run_trna_match.py', config, 'S1', self.write_vcf(d, 'SRC_CHROM=species;SRC_POS=10;SRC_REF=A;SRC_ALT=C'), output)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(record_info(output)['MTTRNA_REGION_MATCH'], 'yes')

    def test_rrna_interval_mode_skips_missing_structure_table_and_enabled_mode_fails_clearly(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); human = d / 'human.tsv'; species = d / 'species.tsv'; missing = d / 'missing.tsv'
            header = 'chrom\tstart\tend\trrna_gene\tstrand\n'
            human.write_text(header + 'chrM\t1\t100\tMT-RNR1\t+\n')
            species.write_text(header + 'species\t1\t100\tMT-RNR1\t+\n')
            input_vcf = self.write_vcf(d, 'SRC_CHROM=species;SRC_POS=10;SRC_REF=A;SRC_ALT=C')
            def config(enabled):
                path = d / f'rrna-{enabled}.yaml'
                path.write_text(f'''rrna_match:
  paths:
    input_vcf_dir: {d}
    fallback_codon_vcf_dir: {d}
    fallback_raw_vcf_dir: {d}
    output_dir: {d}
    reports_dir: {d / 'reports'}
    coordinate_map_dir: {d / 'maps'}
    human_rrna_table: {human}
    species_rrna_table: {species}
  settings:
    input_vcf_pattern: "{{sample}}.vcf"
    fallback_codon_vcf_pattern: "{{sample}}.vcf"
    fallback_raw_vcf_pattern: "{{sample}}.vcf"
    output_suffix: ".out.vcf"
    use_rrna_structure_table: {str(enabled).lower()}
    human_rrna_structure_table: {missing}
    require_pair_pos_match_for_high_conf_stem: true
''')
                return path
            output = d / 'rrna.vcf'
            result = run_script('run_rrna_match.py', config(False), 'S1', input_vcf, output)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(record_info(output)['MTRRNA_STATUS'], 'OK')
            enabled = run_script('run_rrna_match.py', config(True), 'S1', input_vcf, d / 'should-not-exist.vcf')
            self.assertNotEqual(enabled.returncode, 0)
            self.assertIn('rRNA structure annotation is enabled but human structure table is missing', enabled.stderr)


if __name__ == '__main__':
    unittest.main()
